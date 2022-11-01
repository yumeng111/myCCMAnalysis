#include <tuple>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdint.h>
#include <string.h>
#include <argagg/argagg.hpp>

#include "CCMAnalysis/CCMUtils/Utility.h"
#include "CCMAnalysis/CCMUtils/MsgLog.h"
#include "CCMAnalysis/CCMUtils/MakeUniquePatch.h"
#include "CCMAnalysis/CCMFramework/CCMTaskManager.h"

#include "TH1.h"
#include "TApplication.h"

#define RETURNCHECK(VARTYPE, VARNAME, ARG) VARTYPE VARNAME;\
    try {\
        VARNAME = ARG;\
    } catch (ExitStatus & status) {\
        success &= status.status == (int)(EXIT_SUCCESS);\
    }

struct ExitStatus {
    ExitStatus(int status) : status(status) {}
    int status;
};

std::tuple<argagg::parser, argagg::parser_results> parse_arguments(int argc, char ** argv) {
    argagg::parser argparser = {{
        {
            "help", {"-h", "--help"},
                "Print help and exit", 0,
        },
        {
            "output_file", {"-o", "--output-file"},
            "Filename for the output. Only root file output is supported.", 1,
        },
        {
            "input_file", {"-i", "--input-file"},
            "Filename for the input. Both root and binary input files are supported.", 1,
        },
        {
            "input_list", {"-I", "--input-list"},
            "Filename for a list of input files. File should contain input file names\n\t in plain text, separated by line-breaks.", 1,
        },
        {
            "config_file", {"-c", "--config-file"},
            "Filename for the xml configuration file.", 1,
        },
        {
            "log_file", {"-l", "--log-file"},
            "Filename for the log file.", 1,
        },
        {
            "debug_level", {"-d", "--debug-level"},
            "Debugging level for log output.", 1,
        },
        {
            "number_of_events", {"-n", "--number-of-events", "--num-events"},
            "Number of events to process.", 1,
        },
        {
            "remove_first_subrun", {"-s", "--remove-first-subrun", "--remove-subrun"},
            "Remove the first subrun.", 0,
        },
    }};
    std::ostringstream usage;
    usage << argv[0] << std::endl;

    argagg::parser_results args;
    try {
        args = argparser.parse(argc, argv);
    } catch (const std::exception& e) {
        argagg::fmt_ostream fmt(std::cerr);
        fmt << usage.str() << argparser << std::endl
            << "Encountered exception while parsing arguments: " << e.what()
            << std::endl;
        throw ExitStatus(EXIT_FAILURE);
    }

    // Print the usage information when we help is requested
    if(args["help"]) {
        std::cerr << argparser;
        throw ExitStatus(EXIT_SUCCESS);
    }

    return std::make_tuple(argparser, args);
}

std::string get_output(argagg::parser & argparser, argagg::parser_results & args) {
    // Require output file
    if(not args["output_file"]) {
        std::cerr << "--output-file required!" << std::endl;
        throw ExitStatus(EXIT_FAILURE);
    }
    std::string output = args["output_file"].as<std::string>("");
    return output;
}

std::vector<std::string> get_input(argagg::parser & argparser, argagg::parser_results & args) {
    // Require input file
    bool have_input = args["input_file"];
    bool have_input_list = args["input_list"];
    if(have_input and have_input_list) {
        std::cerr << "Only one of --input-file or --input-list can be specified!" << std::endl;
        throw ExitStatus(EXIT_FAILURE);
    }
    if((not have_input) and (not have_input_list)) {
        std::cerr << "Must specify either --input-file or --input-list!" << std::endl;
        throw ExitStatus(EXIT_FAILURE);
    }
    std::vector<std::string> input_files;
    if(have_input_list) {
        std::string input = args["input_list"].as<std::string>("");
        try {
            input_files = Utility::IndirectFileList(input);
        } catch (...) {
            std::cerr << "Issue finding input files!" << std::endl;
            throw ExitStatus(EXIT_FAILURE);
        }
    } else if (have_input) {
        std::string input = args["input_file"].as<std::string>("");
        try {
            input_files = Utility::GetListOfFiles(input);
        } catch (...) {
            std::cerr << "Issue finding input file list!" << std::endl;
            throw ExitStatus(EXIT_FAILURE);
        }
    }
    return input_files;
}

std::string get_config(argagg::parser & argparser, argagg::parser_results & args) {
    // Require config file
    if(not args["config_file"]) {
        std::cerr << "--config-file required!" << std::endl;
        throw ExitStatus(EXIT_FAILURE);
    }
    std::string output = args["config_file"].as<std::string>("");
    return output;
}

std::string get_log(argagg::parser & argparser, argagg::parser_results & args) {
    return args["log_file"].as<std::string>("");
}

int get_debug_level(argagg::parser & argparser, argagg::parser_results & args) {
    return args["debug_level"].as<int>(0);
}

int get_number_of_events(argagg::parser & argparser, argagg::parser_results & args) {
    return args["number_of_events"].as<int>(-1);
}

bool get_remove_first_subrun(argagg::parser & argparser, argagg::parser_results & args) {
    return bool(args["remove_first_subrun"]);
}


//main program
int main (int argc, char** argv) {
try {
    std::tuple<argagg::parser, argagg::parser_results> argarg = parse_arguments(argc, argv);
    argagg::parser & argparser = std::get<0>(argarg);
    argagg::parser_results & args = std::get<1>(argarg);

    // set it so square of the sum of the weights
    // is used to calculate the error on any histograms
    // that are generated
    TH1::SetDefaultSumw2(true);

    // set so ROOT will not manage any memory of ROOT objects
    // that are created
    TH1::AddDirectory(0);

    bool success = true;

    RETURNCHECK(std::vector<std::string>, input_files, get_input(argparser, args))
    RETURNCHECK(std::string, output_file, get_output(argparser, args))
    RETURNCHECK(std::string, config_file, get_config(argparser, args))
    std::string log_file = get_log(argparser, args);
    int debug_level = get_debug_level(argparser, args);
    int number_of_events = get_number_of_events(argparser, args);
    bool remove_first_subrun = get_remove_first_subrun(argparser, args);

    if(not success) {
        std::cerr << argparser;
        throw ExitStatus(EXIT_FAILURE);
    }

    if (remove_first_subrun) {
        MsgInfo("Going to remove the first subrun");
        auto it = input_files.begin();
        while ((it = std::find_if(input_files.begin(),input_files.end(),
                        [](const std::string& s) { return s.find("000000") != std::string::npos;})) != input_files.end()) {
            input_files.erase(it);
        }

        if (input_files.empty()) {
            MsgFatal("No more input files after removing those that match the 0th sub run");
        }
    }

    MsgLog::SetGlobalDebugLevel(debug_level);
    MsgLog::SetPrintRepetitions(false);
    if (!log_file.empty()) {
        MsgInfo(MsgLog::Form("Setting output log file %s",log_file.c_str()));
        MsgLog::SetFileOutput(log_file.c_str());
    }

    std::unique_ptr<CCMTaskManager> taskMan =
        std::make_unique<CCMTaskManager>(config_file, input_files, output_file);

    // check to see if process type is event display
    // if so start the TApplication in order to show
    // the event display
    std::unique_ptr<TApplication> app = nullptr;
    if(taskMan->GetTaskConfig().ProcessType() == "EventDisplay") {
        app = std::unique_ptr<TApplication>(new TApplication("EventDisplay", &argc, argv));
    }

    taskMan->Execute(number_of_events);

    taskMan->Terminate(); // Writes data to output file

    if(app != nullptr) {
        app->Run();
    }

    delete MsgLog::Instance();

    return EXIT_SUCCESS;

} catch(ExitStatus const & status) {
    return status.status;
}
}

