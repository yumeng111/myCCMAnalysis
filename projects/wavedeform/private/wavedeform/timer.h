#ifndef wavedeform_timer_H
#define wavedeform_timer_H

#include <sys/time.h>
#include <sys/resource.h>

#include <icetray/I3Logging.h>

class NNLSTimer {
    std::string name;
    double & sys;
    double & user;
    struct rusage stop, start_time;
    bool fail;
    bool destroyed = false;
    bool print_intermediate;
    bool printed = false;
    bool print_at_all = true;
public:
    void start() {
        fail = (getrusage(RUSAGE_SELF, &start_time) == -1);
        destroyed = false;
    }

    NNLSTimer(std::string name, double& s, double& u, bool print_intermediate=false, bool paa=true) : name(name), sys(s), user(u), print_intermediate(print_intermediate), print_at_all(paa) {
        this->start();
    }

    void print() {
        if(not printed) {
            log_notice("%40s: %9.6fs user %9.8fs system", name.c_str(), user, sys);
            printed = true;
        }
    }

    void end() {
        if (getrusage(RUSAGE_SELF, &stop) != -1 && !fail) {
            double u = 0.0;
            double s = 0.0;

            u += (stop.ru_utime.tv_sec - start_time.ru_utime.tv_sec);
            u += double(stop.ru_utime.tv_usec - start_time.ru_utime.tv_usec) / 1E+06;
            user += u;

            s += (stop.ru_stime.tv_sec - start_time.ru_stime.tv_sec);
            s += double(stop.ru_stime.tv_usec - start_time.ru_stime.tv_usec) / 1E+06;
            sys += s;

            if(print_intermediate)
                if(print_at_all)
                    log_notice("%40s: %9.6fs user %9.8fs system", name.c_str(), u, s);
        }
        destroyed = true;
    }

    ~NNLSTimer() {
        if(destroyed) {
            if(not print_intermediate) {
                if(print_at_all)
                    this->print();
            }
        } else {
            this->end();
            if(not print_intermediate) {
                if(print_at_all)
                    this->print();
            }
        }
    }
};

class DurationTimer {
    double & total;
    struct rusage stop, start;
    double max_time;
    bool fail;
public:
    DurationTimer(double & total, double max_time) : total(total), max_time(max_time) {
        fail = (getrusage(RUSAGE_SELF, &start) == -1);
    }

    double elapsed() {
        fail = (getrusage(RUSAGE_SELF, &stop) == -1);
        double elapsed_time = (stop.ru_utime.tv_sec - start.ru_utime.tv_sec)
            + (stop.ru_stime.tv_sec - start.ru_stime.tv_sec)
            + double(stop.ru_utime.tv_usec - start.ru_utime.tv_usec) / 1E+06
            + double(stop.ru_stime.tv_usec - start.ru_stime.tv_usec) / 1E+06;
        return elapsed_time + total;
    }

    bool timeout() {
        return this->elapsed() > max_time;
    }

    ~DurationTimer() {
        total += this->elapsed();
    }
};

#endif // wavedeform_timer_H
