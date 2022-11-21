import numpy as np
import icecube
from icecube import dataio
from icecube import CCMBinary

def subtract_times(t1, t0):
    # Compute the difference of two unsigned time counters with one overflow bit
    # Cover the trivial case
    if t1 == t0:
        return 0

    # We know the 32nd bit is an overflow bit so we must mask it out
    overflow_0 = t0 & (0x1 << 31)
    overflow_1 = t1 & (0x1 << 31)
    different_overflow = overflow_0 ^ overflow_1
    time_0 = np.int32(np.uint32(t0 & 0x7FFFFFFF))
    time_1 = np.int32(np.uint32(t1 & 0x7FFFFFFF))

    if different_overflow:
        # We know which time is greater
        # But one time has crossed the overflow boundary
        if overflow_0:
            # Convert t1 to a negative number in the frame of t0's zero then subtract
            diff = (time_1 - 0x7FFFFFFE) - time_0
        else:
            # Convert t0 to a negative number in the frame of t1's zero then subtract
            diff = time_1 - (time_0 - 0x7FFFFFFE)
    else:
        # We don't know which time is greater
        # Use the difference with the smallest magnitude

        # Cover the trivial case
        diff_0 = time_1 - time_0

        # Cover the case where the difference crosses the overflow boundary
        if time_0 < time_1:
            # Convert t1 to a negative number in the frame of t0's zero then subtract
            diff_1 = (time_1 - 0x7FFFFFFE) - time_0
        else:
            # Convert t0 to a negative number in the frame of t1's zero then subtract
            diff_1 = time_1 - (time_0 - 0x7FFFFFFE)
        if abs(diff_0) < abs(diff_1):
            diff = diff_0
        else:
            diff = diff_1
    return diff

def empty_mask(frame):
    config = frame["CCMDAQConfig"]
    nboards = len(config.digitizer_boards)
    mask = np.zeros(nboards).astype(bool)
    channel_sizes = np.array(list(frame["CCMTriggerReadout"].triggers[0].channel_sizes))
    # board_times = np.array(list(frame["CCMTriggerReadout"].triggers[0].board_times))
    for i in range(nboards):
        mask[i] = np.any(channel_sizes[i*16:(i+1)*16] > 0)
    return mask

class TimeReader:
    def __init__(self, fname):
        self.f = dataio.I3File(fname)
        self.pop()
        time_read = np.array(list(self.frame["CCMTriggerReadout"].triggers[0].board_times))
        computer_time_read = np.array(list(self.frame["CCMTriggerReadout"].triggers[0].board_computer_times))
        mask = empty_mask(self.frame)
        self.time_cache = np.zeros((len(time_read), 0)).tolist()
        self.computer_time_cache = np.zeros((len(time_read), 0)).tolist()
        for i in range(len(self.time_cache)):
            if mask[i]:
                self.time_cache[i].append(time_read[i])
        self.last_raw_times = np.zeros(len(self.time_cache)).astype(np.uint32)
        self.last_times = np.zeros(len(self.time_cache)).astype(np.int64)

    def pop(self):
        if not self.f.more():
            return False
        frame = self.f.pop_daq()
        while "CCMTriggerReadout" not in frame.keys() or len(frame["CCMTriggerReadout"].triggers) == 0:
            if not self.f.more():
                return False
            frame = self.f.pop_daq()
        self.frame = frame
        return True

    def pop_times(self):
        self.pop()
        time_read = np.array(list(self.frame["CCMTriggerReadout"].triggers[0].board_times))
        mask = empty_mask(self.frame)
        for i in range(len(self.time_cache)):
            if mask[i]:
                self.time_cache[i].append(time_read[i])

    def get_times(self):
        times_to_return = np.zeros(len(self.time_cache)).astype(np.int64)
        while not np.all([len(l) > 0 for l in self.time_cache]):
            self.pop_times()
        for i in range(len(times_to_return)):
            raw_time = self.time_cache[i].pop(0)
            time_diff = subtract_times(raw_time, self.last_raw_times[i])
            abs_time = self.last_times[i] + time_diff
            times_to_return[i] = abs_time

            self.last_raw_times[i] = raw_time
            self.last_times[i] = abs_time
        return times_to_return

def find_pairs(times0, times1, offset1, max_delta):
    good_pairs = []
    orphans = []
    pairs = []
    i = 0
    j = 0
    last_i = None
    last_j = None
    while i < len(times0) or j < len(times1):
        if i == len(times0):
            orphaned_pair = (None, j)
            pairs.append(orphaned_pair)
            orphans.append(orphaned_pair)
            j += 1
            continue
        elif j == len(times1):
            orphaned_pair = (i, None)
            pairs.append(orphaned_pair)
            orphans.append(orphaned_pair)
            i += 1
            continue
        time_diff = times0[i] - (times1[i] + offset1)
        if abs(time_diff) <= max_delta:
            if (last_i is not None and last_i < i-1) or (last_j is not None and last_j < j-1):
                triplets = []
                for _i in range(last_i + 1, i):
                    triplets.append((_i, None, times0[_i]))
                for _j in range(last_j + 1, j):
                    triplets.append((None, _j, times1[_j] + offset1))
                triplets = sorted(triplets, key=lambda x: x[2])
                orphaned_pairs = [t[:2] for t in triplets]
                pairs.extend(orphaned_pairs)
                orphans.extend(orphaned_pairs)

            pairs.append((i, j))
            good_pairs.append(pairs[-1])
            last_i = i; last_j = j
            i += 1; j += 1
        else:
            if times0[i] > times1[j] + offset1:
                i += 1
            else:
                j += 1

    return pairs, good_pairs, orphans

class TimeComparison:
    def __init__(self, fname0, fname1):
        self.reader_0 = TimeReader(fname0)
        self.reader_1 = TimeReader(fname1)

m_reader = TimeReader("mills.i3.zst")
m_times = np.array([m_reader.get_times() for i in range(100)]).T
w_reader = TimeReader("wills.i3.zst")
w_times = np.array([w_reader.get_times() for i in range(100)]).T

for i in range(1, len(w_times)):
    pairs, good_pairs, orphans = find_pairs(w_times[0], w_times[i], w_times[0][0] - w_times[i][0], 16)
    print(len(pairs), len(good_pairs), len(orphans))

#for delta_trigger in range(-50, 51):
#    find_pairs(m_times[0], w_times[1])


