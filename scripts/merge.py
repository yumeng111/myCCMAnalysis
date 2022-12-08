import numpy as np
import icecube
import I3Tray
from icecube import icetray
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
    last_idx = 0
    for i in range(nboards):
        n_channels = len(config.digitizer_boards[i].channels)
        next_idx = last_idx + n_channels
        mask[i] = np.any(channel_sizes[last_idx:next_idx] > 0)
        last_idx = next_idx
    return mask

class TimeReader:
    def __init__(self, fnames, n_skip=0):
        self.f = dataio.I3FrameSequence(fnames)
        self.pop()
        time_read = np.array(list(self.frame["CCMTriggerReadout"].triggers[0].board_times))
        computer_time_read = np.array(list(self.frame["CCMTriggerReadout"].triggers[0].board_computer_times))
        mask = empty_mask(self.frame)
        self.time_cache = np.zeros((len(time_read), 0)).tolist()
        self.last_raw_time = np.zeros(len(time_read)).astype(np.uint32)
        for i in range(len(self.time_cache)):
            if mask[i]:
                self.time_cache[i].append(time_read[i])
                self.last_raw_time[i] = time_read[i]
        while not np.all([len(self.time_cache[i]) > n_skip for i in range(len(self.time_cache))]):
            self.pop_times()
        for i in range(len(self.time_cache)):
            if len(self.time_cache[i]) > 0:
                self.time_cache[i] = self.time_cache[i][n_skip:]

    def n_boards(self):
        return len(self.time_cache)

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
                raw_time = time_read[i]
                abs_time = self.time_cache[i][-1] + subtract_times(raw_time, self.last_raw_time[i])
                self.time_cache[i].append(abs_time)
                self.last_raw_time[i] = raw_time

    def get_times(self, N, board_idx):
        while len(self.time_cache[board_idx]) < N:
            self.pop_times()
        return list(self.time_cache[board_idx][:N])

def find_pairs(times0, times1, offset1, max_delta):
    good_pairs = []
    orphans = []
    pairs = []
    i = 0
    j = 0
    last_i = None
    last_j = None
    while i < len(times0) or j < len(times1):
        if i == len(times0) or j == len(times1):
            i = len(times0)
            j = len(times1)
            if last_i is None:
                last_i = -1
            if last_j is None:
                last_j = -1
            triplets = []
            for _i in range(last_i + 1, i):
                triplets.append((_i, None, times0[_i]))
            for _j in range(last_j + 1, j):
                triplets.append((None, _j, times1[_j] + offset1))
            triplets = sorted(triplets, key=lambda x: x[2])
            orphaned_pairs = [t[:2] for t in triplets]
            pairs.extend(orphaned_pairs)
            orphans.extend(orphaned_pairs)
            last_j = j
            last_i = i
            continue

        time_diff = times0[i] - (times1[j] + offset1)
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
            if times0[i] < times1[j] + offset1:
                i += 1
            else:
                j += 1

    return pairs, good_pairs, orphans

def get_time_delta(times0, times1, delta_trigger):
    if delta_trigger < 0:
        delta = times0[0] - times1[-delta_trigger]
    else:
        delta = times0[delta_trigger] - times1[0]
    return delta

def compute_trigger_offset(reader0, board_idx0, reader1, board_idx1, jitter_tests=None, max_delta=2, min_triggers=100, max_triggers=2000, increment=25, threshold=0.9):
    if jitter_tests is None:
        jitter_tests = [-2, 2]
    n_triggers = min_triggers
    prev_min = 0
    prev_max = 0
    new_min = -int(n_triggers/2)
    new_max = int(n_triggers/2)+1
    times0 = np.array(reader0.get_times(n_triggers, board_idx0))
    times1 = np.array(reader1.get_times(n_triggers, board_idx1))
    deltas_above = list(range(prev_max, new_max))
    deltas_below = list(range(new_min, prev_min))
    prev_min = new_min
    prev_max = new_max
    results = []

    while True:
        for delta_trigger in deltas_above + deltas_below:
            delta = get_time_delta(times0, times1, delta_trigger)
            pairs, good_pairs, orphans = find_pairs(times0, times1, delta, max_delta)
            results.append((delta_trigger, delta, len(pairs), len(good_pairs), len(orphans)))
            for t in jitter_tests:
                pairs, good_pairs, orphans = find_pairs(times0, times1, delta + t, max_delta)
                results.append((delta_trigger, delta + t, len(pairs), len(good_pairs), len(orphans)))

        best_pair_result = max(results, key=lambda x: x[3])
        delta_trigger, delta, pairs, good_pairs, orphans = best_pair_result

        print(len(times0))
        print("With", n_triggers, "triggers found", good_pairs, "good pairs")
        print("With an offset of", delta_trigger, "we expect at most", (n_triggers - abs(delta_trigger)), "good pairs")
        print("Achieved a ratio of", good_pairs / (n_triggers - abs(delta_trigger)))
        print()
        if good_pairs / (n_triggers - abs(delta_trigger)) >= threshold:
            break

        n_triggers += increment
        if(n_triggers > max_triggers):
            return None
        new_max = int(0.5*n_triggers)+1
        new_min = -int(0.5*n_triggers)
        deltas_above = list(range(prev_max, new_max))
        deltas_below = list(range(new_min, prev_min))
        times0 = np.array(reader0.get_times(n_triggers, board_idx0))
        times1 = np.array(reader1.get_times(n_triggers, board_idx1))
        prev_min = new_min
        prev_max = new_max

    return best_pair_result, n_triggers

class MergedSource(icetray.I3Module) :
    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddParameter("FileLists", "File lists to merge", [])
        self.AddParameter("MaxTimeDiff", "Maximum time difference between associated triggers (ns)", 16)

    def Configure(self):
        #print("Configure")
        file_lists = self.GetParameter("FileLists")
        max_delta_t = self.GetParameter("MaxTimeDiff")
        self.max_delta = int(np.ceil(max_delta_t / 8))
        n_daqs = len(file_lists)
        readers = [TimeReader(file_lists[i]) for i in range(n_daqs)]
        self.n_boards = [r.n_boards() for r in readers]
        offsets = [[None for i in range(n)] for n in self.n_boards]
        offsets[0][0] = 0
        for i in range(1, len(offsets)):
            result = compute_trigger_offset(readers[0], 0, readers[i], 0)
            if result is None:
                s = "Error: cannot align DAQ 0 and DAQ " + str(i)
                #print(s)
                raise RuntimeError(s)
            else:
                (delta_trigger, delta, pairs, good_pairs, orphans), n_triggers = result
                offsets[i][0] = delta
        for i in range(len(offsets)):
            for j in range(1, self.n_boards[i]):
                result = compute_trigger_offset(readers[0], 0, readers[i], j)
                (delta_trigger, delta, pairs, good_pairs, orphans), n_triggers = result
                offsets[i][j] = delta
        del readers

        self.frame_sequences = [dataio.I3FrameSequence(file_lists[i]) for i in range(n_daqs)]

        self.offsets = offsets
        self.frame_cache = np.zeros((len(self.offsets), 0)).tolist()
        self.frame_idxs = [np.zeros(n).astype(int)-1 for n in self.n_boards]
        self.empty_cache = [np.zeros(n).astype(bool) for n in self.n_boards]
        self.mask_cache = np.zeros((len(self.offsets), 0)).tolist()
        self.time_cache = [np.zeros((n, 0)).tolist() for n in self.n_boards]
        self.last_raw_time = [np.zeros(n).astype(np.uint32) for n in self.n_boards]
        self.last_time = [np.zeros(n).astype(np.uint64) for n in self.n_boards]
        self.configs = [None for n in self.n_boards]
        self.fill_computer_time = None
        self.push_config = False
        self.frames_to_push = []
        self.next_triggers()

    def pop_frame(self, idx):
        #print("pop_frame")
        if not self.frame_sequences[idx].more():
            return False
        frame = self.frame_sequences[idx].pop_frame()
        while frame.Stop != icecube.icetray.I3Frame.DAQ:
            if frame.Stop == icecube.icetray.I3Frame.Geometry:
                self.configs[idx] = frame["CCMDAQConfig"]
                self.push_config = True
            elif frame.Stop != icecube.icetray.I3Frame.DAQ:
                self.frames_to_push.append(frame)
            if not self.frame_sequences[idx].more():
                return False
            frame = self.frame_sequences[idx].pop_frame()
        while "CCMTriggerReadout" not in frame.keys() or len(frame["CCMTriggerReadout"].triggers) == 0:
            if not self.frame_sequences[idx].more():
                return False
            frame = self.frame_sequences[idx].pop_frame()
            while frame.Stop != icecube.icetray.I3Frame.DAQ:
                if frame.Stop == icecube.icetray.I3Frame.Geometry:
                    self.configs[idx] = frame["CCMDAQConfig"]
                    self.push_config = True
                elif frame.Stop != icecube.icetray.I3Frame.DAQ:
                    self.frames_to_push.append(frame)
                if not self.frame_sequences[idx].more():
                    return False
                frame = self.frame_sequences[idx].pop_frame()
        self.frame_cache[idx].append(frame)
        mask = empty_mask(frame)
        self.mask_cache[idx].append(mask)
        time_read = np.array(list(frame["CCMTriggerReadout"].triggers[0].board_times))
        if self.fill_computer_time is None:
            self.fill_computer_time = len(frame["CCMTriggerReadout"].triggers[0].board_computer_times) > 0
        for i in range(len(self.time_cache[idx])):
            if mask[i]:
                raw_time = time_read[i]
                rel_time = subtract_times(raw_time, self.last_raw_time[idx][i])
                if len(self.time_cache[idx][i]):
                    abs_time = self.last_time[idx][i] + rel_time
                else:
                    abs_time = rel_time
                self.time_cache[idx][i].append(abs_time)
                self.last_raw_time[idx][i] = raw_time
                self.last_time[idx][i] = abs_time
            else:
                self.time_cache[idx][i].append(np.inf)
        return True

    def get_config_frame(self):
        config = CCMBinary.CCMDAQConfig()
        for c in self.configs:
            config.machine_configurations.extend(c.machine_configurations)
            config.digitizer_boards.extend(c.digitizer_boards)
        frame = icecube.icetray.I3Frame(icecube.icetray.I3Frame.Geometry)
        frame.Put("CCMDAQConfig", config, icecube.icetray.I3Frame.Geometry)
        return frame

    def next_trigger(self, daq_idx, board_idx):
        #print()
        #print("next_trigger")
        current_frame_idx = self.frame_idxs[daq_idx][board_idx]
        #print("frame_cache_length", len(self.frame_cache[daq_idx]))
        #print("current_frame_idx", current_frame_idx)
        frame_idx = current_frame_idx + 1
        #print("target_frame_idx", frame_idx)
        while True:
            while len(self.frame_cache[daq_idx]) <= frame_idx:
                res = self.pop_frame(daq_idx)
                if not res:
                    return False
            if self.mask_cache[daq_idx][frame_idx][board_idx]:
                break
            frame_idx += 1
        #print("result_frame_idx", frame_idx)
        #print("frame_cache_length", len(self.frame_cache[daq_idx]))
        #print()
        self.frame_idxs[daq_idx][board_idx] = frame_idx
        return True

    def next_triggers(self):
        #print("next_triggers")
        for daq_idx in range(len(self.frame_idxs)):
            for board_idx in range(self.n_boards[daq_idx]):
                res = self.next_trigger(daq_idx, board_idx)
                if not res:
                    return False
        self.clear_unused_frames()
        return True

    def clear_unused_frames(self):
        #print()
        #print("clear_unused_frames")
        #print(self.frame_idxs)
        #print([[len(x) for x in l] for l in self.time_cache])
        #print([len(l) for l in self.frame_cache])
        for daq_idx in range(len(self.frame_cache)):
            min_idx = np.amin(self.frame_idxs[daq_idx])
            for i in range(min_idx):
                self.frame_cache[daq_idx].pop(0)
                self.mask_cache[daq_idx].pop(0)
                for j in range(self.n_boards[daq_idx]):
                    self.time_cache[daq_idx][j].pop(0)
            self.frame_idxs[daq_idx] -= min_idx
        #print([[len(x) for x in l] for l in self.time_cache])
        #print([len(l) for l in self.frame_cache])
        #print()

    def get_trigger_readout(self):
        #print("get_trigger_readout")
        readout = CCMBinary.CCMTriggerReadout()
        readout.triggers.append(CCMBinary.CCMTrigger())

        #print(self.frame_idxs)
        #print([[len(x) for x in l] for l in self.time_cache])
        times = [[self.time_cache[daq_idx][board_idx][self.frame_idxs[daq_idx][board_idx]] + self.offsets[daq_idx][board_idx] for board_idx in range(n)] for daq_idx, n in enumerate(self.n_boards)]
        min_time = min([np.amin(times[i]) for i in range(len(self.n_boards))])

        all_bad = np.all([np.all(self.empty_cache[daq_idx]) for daq_idx in range(len(self.n_boards))])

        if all_bad or min_time == np.inf:
            return None

        for daq_idx, n in enumerate(self.n_boards):
            last_idx = 0
            for board_idx in range(n):
                n_channels = len(self.configs[daq_idx].digitizer_boards[board_idx].channels)
                next_idx = last_idx + n_channels
                if times[daq_idx][board_idx] - min_time > self.max_delta or self.empty_cache[daq_idx][board_idx]:
                    CCMBinary.merge_empty_trigger(readout, n_channels, self.fill_computer_time)
                else:
                    tr = self.frame_cache[daq_idx][self.frame_idxs[daq_idx][board_idx]]["CCMTriggerReadout"]
                    CCMBinary.merge_triggers(readout, tr, board_idx, last_idx, next_idx, self.fill_computer_time)
                    res = self.next_trigger(daq_idx, board_idx)
                    if not res:
                        times[daq_idx][board_idx] = np.inf
                        self.empty_cache[daq_idx][board_idx] = True
                last_idx = next_idx
        self.clear_unused_frames()
        return readout

    def Process(self):
        if self.push_config:
            self.PushFrame(self.get_config_frame())
            self.push_config = False
        frame = icecube.icetray.I3Frame(icecube.icetray.I3Frame.DAQ)
        readout = self.get_trigger_readout()
        if readout is None:
            self.RequestSuspension()
            return
        frame.Put("CCMDigitalReadout", CCMBinary.I3VectorI3VectorUInt16(readout.samples), icecube.icetray.I3Frame.DAQ)
        frame.Put("CCMTriggers", CCMBinary.I3VectorCCMTrigger(readout.triggers), icecube.icetray.I3Frame.DAQ)
        self.PushFrame(frame)

reader0 = TimeReader(["mills.i3.zst"])
reader1 = TimeReader(["wills.i3.zst"])
best_pair_result, n_triggers = compute_trigger_offset(reader0, 0, reader1, 0)
delta_trigger, delta, pairs, good_pairs, orphans = best_pair_result

print("Found best match with", best_pair_result[2], "good pairs and", best_pair_result[3], "orphans after testing", n_triggers, "triggers")
print("Offset by", ("+"+str(delta_trigger)) if delta_trigger > 0 else delta_trigger, "triggers")
print("With a time offset of", ("+"+str(-delta)) if -delta > 0 else -delta, "time steps")

tray = I3Tray.I3Tray()

tray.AddModule(MergedSource, "reader", FileLists=[["mills.i3.zst"],["wills.i3.zst"]])
tray.AddModule("I3MultiWriter", "writer", SizeLimit=1, SyncStream=icetray.I3Frame.Geometry, FileName="test-%04d.i3.zst", MetadataStreams=[icetray.I3Frame.Geometry, icetray.I3Frame.Calibration, icetray.I3Frame.DetectorStatus])

tray.Execute()

# reader0 = TimeReader("mills.i3.zst")
# reader1 = TimeReader("wills.i3.zst")
# times0 = reader0.get_times()
# times1 = reader1.get_times()
# real_time_diff = times0[0] - times1[0]
# print("Real time diff", real_time_diff)
# reader0 = TimeReader("mills.i3.zst")
# reader1 = TimeReader("wills.i3.zst")
# for i in range(10):
#     reader1.get_times()
# min_triggers = 100
# times0 = np.array([reader0.get_times() for i in range(min_triggers)]).T
# times1 = np.array([reader1.get_times() for i in range(min_triggers)]).T
# delta_trigger = 10
# delta = get_time_delta(times0, times1, delta_trigger)
# print("Test time diff", delta)
# pairs, good_pairs, orphans = find_pairs(times0[0], times1[0], delta, 2)
# print(len(pairs), len(good_pairs), len(orphans))

# for delta_trigger in range(-50, 51):
#     delta = get_time_delta(times0, times1, delta_trigger)
#     pairs, good_pairs, orphans = find_pairs(times0[0], times1[0], delta, 2)
#     print(delta_trigger, len(pairs), len(good_pairs), len(orphans))
