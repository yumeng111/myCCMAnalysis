#ifndef MSGLOG_H
#define MSGLOG_H
/*------------------------------------------------------
 *   MsgLog
 *    Author: J.L. Klay (CalPoly)                           
 *    Date: 29-Apr-2008
 *-------------------------------------------------------*/

#include <string>
#include <iostream>

class MsgLog {

  public:
    MsgLog();  //default constructor
    ~MsgLog();  //destructor

    static MsgLog* Instance() { return fgInstance; }

    enum EType {kFatal = 0, kError, kWarning, kInfo, kDebug, kMaxType};
    typedef void (*MsgLogNotification)(EType type, const char* message);

    static void  EnableDebug(bool enabled);
    static void  SetGlobalLogLevel(EType type);
    static int   GetGlobalLogLevel();
    static void  SetGlobalDebugLevel(int level);
    static int   GetGlobalDebugLevel();

    static void  SetStandardOutput();
    static void  SetStandardOutput(EType type);
    static void  SetErrorOutput();
    static void  SetErrorOutput(EType type);
    static void  SetFileOutput(const char* fileName);
    static void  SetFileOutput(EType type, const char* fileName);
    static void  SetStreamOutput(std::ostream* stream);
    static void  SetStreamOutput(EType type, std::ostream* stream);
    static void  SetLogNotification(MsgLogNotification pCallBack);
    static void  SetLogNotification(EType type, MsgLogNotification pCallBack);
    static void  Flush();

    static void  SetPrintType(bool on);
    static void  SetPrintType(EType type, bool on);
    static void  SetPrintLocation(bool on);
    static void  SetPrintLocation(EType type, bool on);

    static void  SetPrintRepetitions(bool on);

    static const std::string Form(const char *fmt, ...);

    // the following public methods are used by the preprocessor macros
    // and should not be called directly
    static bool IsDebugEnabled() {return fgDebugEnabled;}
    static void Message(int level, const std::string message,
        const char* function, const char* file, int line);
    static void Debug(int level, const std::string message,
        const char* function, const char* file, int line);

    static int  RedirectStdoutTo(EType type, int level, const char* function,
        const char* file, int line, bool print);
    static int  RedirectStderrTo(EType type, int level, const char* function,
        const char* file, int line, bool print);
    static void RestoreStdout(int original);
    static void RestoreStderr(int original);

    static std::ostream& Stream(EType type, int level,
        const char* function, const char* file, int line);

  private:
    MsgLog(const MsgLog& log); //copy constructor
    MsgLog& operator = (const MsgLog& log);

    void ReadEnvSettings();

    void CloseFile(int type);
    FILE* GetOutputStream(int type);

    void  PrintMessage(int type, const std::string message,
        const char* function,
        const char* file, int line);

    void  PrintString(int type, FILE* stream, const char* format, ...);
    void  PrintRepetitions();

    const char* Test(const char* file, const char* function);

    int   RedirectTo(FILE* stream, EType type, int level,
        const char* function,
        const char* file, int line, bool print);

    std::ostream&  GetStream(EType type, int level,
        const char* function, const char* file, int line);


    enum {kDebugOffset = kDebug-1};

    static MsgLog* fgInstance;      //pointer to current instance

    static bool       fgDebugEnabled;  //flag for debug en-/disabling

    int               fGlobalLogLevel; // global logging level

    int               fOutputTypes[kMaxType];     // types of output streams
    std::string       fFileNames[kMaxType];       // file names
    FILE*             fOutputFiles[kMaxType];     // log output files
    std::ostream*     fOutputStreams[kMaxType];   // log output streams

    bool              fPrintType[kMaxType];       // print type on/off
    bool              fPrintLocation[kMaxType];   // print file and line on/off

    bool              fPrintRepetitions;          // print number of repetitions instead of
    // repeated message on/off

    int               fRepetitions;               // counter of repetitions
    int               fLastType;                  // type of last message
    std::string       fLastMessage;               // last message
    std::string       fLastFunction;              // function name of last message
    std::string       fLastFile;                  // file name of last message
    int               fLastLine;                  // line number of last message
    MsgLogNotification fCallBacks[kMaxType];   // external notification callback

};

// function name
#ifdef __FUNCTION__
#define FUNCTIONNAME() __FUNCTION__
#else
#define FUNCTIONNAME() __func__
#endif

// redirection
#define REDIRECTSTDOUT(type, level, whatever) {int originalStdout = MsgLog::RedirectStdoutTo(type, level, FUNCTIONNAME(), __FILE__, __LINE__, kFALSE); whatever; MsgLog::RestoreStdout(originalStdout);}
#define REDIRECTSTDERR(type, level, whatever) {int originalStderr = MsgLog::RedirectStderrTo(type, level, FUNCTIONNAME(), __FILE__, __LINE__, kFALSE); whatever; AliLog::RestoreStderr(originalStderr);}
#define REDIRECTSTDOUTANDSTDERR(type, level, whatever) {int originalStdout = MsgLog::RedirectStdoutTo(type, level, FUNCTIONNAME(), __FILE__, __LINE__, kFALSE); int originalStderr = MsgLog::RedirectStderrTo(type, level, FUNCTIONNAME(), __FILE__, __LINE__, kFALSE); whatever; MsgLog::RestoreStderr(originalStderr); MsgLog::RestoreStdout(originalStdout);}

// debug level
#ifdef LOG_NO_DEBUG
#define MsgDebugLevel() -1
#else
#define MsgDebugLevel() (MsgLog::IsDebugEnabled())
#endif

// debug messages
#ifdef LOG_NO_DEBUG
#define MsgDebug(level, message)
#else
#define MsgDebug(N, A) \
  do { \
    if (!MsgLog::IsDebugEnabled()) break; \
    MsgLog::Debug(N, A, FUNCTIONNAME(), \
        __FILE__, __LINE__); } while (false)
#endif

// redirection to debug
#define StdoutToMsgDebug(level, whatever) REDIRECTSTDOUT(MsgLog::kDebug, level, whatever)
#define StderrToMsgDebug(level, whatever) REDIRECTSTDERR(MsgLog::kDebug, level, whatever)
#define ToMsgDebug(level, whatever) REDIRECTSTDOUTANDSTDERR(MsgLog::kDebug, level, whatever)

//debug stream objects
#define MsgDebugStream(level) MsgLog::Stream(MsgLog::kDebug, level, FUNCTIONNAME(), __FILE__, __LINE__)

// info messages
#ifdef LOG_NO_INFO
#define MsgInfo(message)
#else
#define MsgInfo(message) {MsgLog::Message(MsgLog::kInfo, message, FUNCTIONNAME(), __FILE__, __LINE__);}
#endif

//redirection to info
#define StdoutToMsgInfo(whatever) REDIRECTSTDOUT(MsgLog::kInfo, 0, whatever)
#define StderrToMsgInfo(whatever) REDIRECTSTDERR(MsgLog::kInfo, 0, whatever)
#define ToMsgInfo(whatever) REDIRECTSTDOUTANDSTDERR(MsgLog::kInfo, 0, whatever)

//info stream objects
#define MsgInfoStream() MsgLog::Stream(MsgLog::kInfo, 0, FUNCTIONNAME(), __FILE__, __LINE__)

// warning messages
#ifdef LOG_NO_WARNING
#define MsgWarning(message)
#else
#define MsgWarning(message) {MsgLog::Message(MsgLog::kWarning, message, FUNCTIONNAME(), __FILE__, __LINE__);}
#endif

// redirection to warning
#define StdoutToMsgWarning(whatever) REDIRECTSTDOUT(MsgLog::kWarning, 0, whatever)
#define StderrToMsgWarning(whatever) REDIRECTSTDERR(MsgLog::kWarning, 0, whatever)
#define ToMsgWarning(whatever) REDIRECTSTDOUTANDSTDERR(MsgLog::kWarning, 0, whatever)

// warning stream objects
#define MsgWarningStream() MsgLog::Stream(MsgLog::kWarning, 0, FUNCTIONNAME(), __FILE__, __LINE__)

// error messages
#define MsgError(message) {MsgLog::Message(MsgLog::kError, message, FUNCTIONNAME(), __FILE__, __LINE__);}

// redirection to error
#define StdoutToMsgError(whatever) REDIRECTSTDOUT(MsgLog::kError, 0, whatever)
#define StderrToMsgError(whatever) REDIRECTSTDERR(MsgLog::kError, 0, whatever)
#define ToMsgError(whatever) REDIRECTSTDOUTANDSTDERR(MsgLog::kError, 0, whatever)

// error stream objects
#define MsgErrorStream() MsgLog::Stream(MsgLog::kError, 0, FUNCTIONNAME(), __FILE__, __LINE__)

// fatal messages                                                                                                                       
#define MsgFatal(message) {MsgLog::Message(MsgLog::kFatal, message, FUNCTIONNAME(), __FILE__, __LINE__);}

#endif //MSGLOG_H
