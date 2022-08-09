#ifndef CCMROOTIO_H
#define CCMROOTIO_H
/*!
 *
 * \file CCMRootIO.h
 * \brief Handles the input/output of ROOT objects.
 * \note Adapted from SBRootIO of SciBath,
 * \author R. T. Thornton (LANL)
 * \date May 7, 2020
 */

#include <memory>
#include <random>
#include <vector>
#include <cstring>
#include <stdint.h>

#include "CCMAnalysis/CCMIO/CCMEventTreeHandle.h"
#include "CCMAnalysis/CCMUtils/Utility.h"

class TFile;
class RawData;
class Events;
class AccumWaveform;
class Pulses;
class MCTruth;

/*!
 * \class CCMRootIO
 * \brief Handles the input/output of ROOT objects.
 * \note Adapted from SBRootIO of SciBath,
 * \author R. T. Thornton (LANL)
 * \date May 7, 2020
 */
class CCMRootIO 
{
  public:
    /*!
     * \fn CCMRootIO
     * \brief constructor
     */
    CCMRootIO();

    /*!
     * \fn virtual ~CCMRootIO
     * \brief virtual destructor
     */
    virtual ~CCMRootIO();

    //Methods related to the input file list
    /*!
     * \fn int AddFile(const char* file_regexp, bool hasWildcards=false)
     * \brief Adds file_regexp to the list of input files
     * \param[in] file_regexp name of the file
     * \param[in] hasWildcards boolean to see if extra info is given in title 
     * \return The total number of files in list
     */
    int AddFile(const char* file_regexp, bool hasWildcards=false);
    /*!
     * \fn int RemoveFile(const char* file_regexp)
     * \brief Removes file from input list matching file_regexp name
     * \param[in] file_regexp file to remove from list
     * \return An integer
     * \warning File not implemented yet
     */
    int RemoveFile(const char* file_regexp);
    /*!
     * \fn int GoToFile(const char* file)
     * \brief Go to file with name matching file
     * \param[in] file the name of the file to go to
     * \return 1 if found and 0 if not found
     */
    int GoToFile(const char* file);
    /*!
     * \fn int AdvanceFile(int n=1)
     * \brief Advance to the nth file from the current file
     * \param[in] n the number of files to advance
     * \return The number of files advanced
     */
    int AdvanceFile(int n=1);
    /*!
     * \fn int RewindFile(int n=1)
     * \brief Rewind to the nth file from the current file
     * \param[in] n the number of files to skip
     * \return The number of files skipped
     */
    int RewindFile(int n=1);
    /*!
     * \fn const char* CurrentFileName() const
     * \brief Returns the name of the current file as a const char pointer
     */
    const char* CurrentFileName()   const;
    /*!
     * \fn const char* FileName(int i) const
     * \brief Returns the name of the ith file as a const char pointer
     * \param[in] i the index of which file to look at
     */
    const char* FileName(int i) const;

    /*!
     * \fn void SetInFileName(const char* infile)
     * \brief When only using one file, one can set the name by using this function
     * \param[in] infile name of the infile
     */
    void SetInFileName(std::string const & infile);
    /*!
     * \fn SetInFileList(std::vector<std::string> infileList)
     * \brief Set the entire file list at once
     * \param[in] infileList a vector of strings containing the names of all the input files
     */
    void SetInFileList(std::vector<std::string> infileList);
    /*!
     * \fn int SetupInputFile()
     * \brief Sets up the file for ussage
     * 
     * Declares the TFile and makes tells all the trees to 
     * do the initial grab of all the branches
     *
     * \return nonzero of successful 
     */
    int SetupInputFile();

    // Methods related to the output file
    /*!
     * \fn void SetOutFileName(std::string outfile)
     * \brief Sets the name of the outfile
     * \param[in] outfile name of the outfile
     */
    void SetOutFileName(std::string outfile) { fOutFileName = outfile; }
    /*!
     * \fn void SetupOutputFile()
     * \brief Same as SetupInputFile(), but for the output file
     */
    void SetupOutputFile();
    /*!
     * \fn const char* OutFileName() const
     * \brief Returns the name of the output file as a const char pointer
     */
    const char* OutFileName() const { return fOutFileName.c_str(); }

    /*!
     * \fn void Close()
     * \brief Closes all files and writes anything unwriten to the output file
     */
    void Close();

    /*!
     * \fn void Clear()
     * \brief Clears the input file list
     */
    void Clear();

    /*!
     * \fn bool ReadOK()
     * \brief Returns true if the last read from the event tree was sucessful
     */
    bool ReadOK() const { return fReadOK; }

    //void SetOutSizeLimit(int mbLimit) {
    //  fOutSizeLimit = mbLimit; fOutSizeLimit *= 1000000;
    //}

    /*!
     * \fn CCMEventTreeHandle& GetEventTree()
     * \brief Returns current reference variable to the event tree handle class
     */
    CCMEventTreeHandle& GetEventTree();

    bool ReadTree() { return fEventHandle->ReadTree(); }

    AccumWaveform& GetAccumWaveform();
    Events& GetEvents();
    RawData& GetRawData();
    Pulses& GetPulses();
    MCTruth& GetMCTruth();

    void SetAccumWaveform(const AccumWaveform & accumWaveform);
    void SetEvents(const Events & event);
    void SetRawData(const RawData & rawData);
    void SetPulses(const Pulses & pulses);
    void SetMCTruth(const MCTruth & mcTruth);

    void Dump();

    // Methods related to input event list
    uint32_t GoTo(uint32_t event);
    uint32_t Advance(uint32_t n = 1);
    uint32_t Rewind(uint32_t n = 1);
    int Reload();
    int WriteTrigger();

    uint32_t NumInputFiles() { return fInFileList.size(); }
    uint32_t GetNumOfEvents(std::string fileName = "");
    uint32_t GoToRandom();

    uint32_t GetTriggerNumber() const { return fTriggerNumber; }

    void SetParameter(std::string name, const int value);
    void SetParameter(std::string name, const double value);
    void SetParameter(std::string name, std::string value);


  protected:
    virtual void UpdateTriggerNumbers();

    bool fReadOK;         ///< Next read should be OK?
    uint32_t fTriggerNumber;    ///< Trigger number for current event

    std::unique_ptr<CCMEventTreeHandle> fEventHandle;     ///< The break in board tree handle

  private:

    //Input stream data
    int fFileIndex;              ///< Current place in the file list
    std::shared_ptr<TFile> fInFile;                 ///< Input file pointer
    bool fOwnHandle;
    std::vector<std::string> fInFileList;       ///< List of files attached
    std::vector<uint32_t> fInFileEntries;
    std::vector<uint32_t> fInFileEntriesCDF;

    //Output stream data
    std::string fOutFileName;                   ///< Name of the output file
    std::shared_ptr<TFile> fOutFile;                       ///< Output data file pointer
    uint32_t    fNWrite;                ///< Number of events written
    uint32_t    fFlushFreq;                     ///< Flush output every n events
    long long   fOutSizeLimit;                  ///< Output size limit

    std::random_device fRD;
    std::mt19937_64 fMT;
};

#endif // CCMROOTIO_H
