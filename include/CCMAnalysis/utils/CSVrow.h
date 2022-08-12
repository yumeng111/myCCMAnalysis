/*!**********************************************
 * \file CSVrow.h
 * \author R. T. Thornton (LANL)
 * \date February 24, 2020
 * \brief Contains CSVRow and CSVIterator classes to parse a .csv file
 *
 * Contains CSVRow and CSVIterator classes to parse  and store a .csv file. 
 * copied from https://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c
 *
 ***********************************************/
#ifndef CSVrow_h
#define CSVrow_h

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>

//-------------------------------------------------------------------------------------------------
/*!**********************************************
 * \class CSVRow
 * \brief Container to store the contents of the csv file
 ***********************************************/
class CSVRow
{
  public:
    /// \fn std::string const & operator[] (std::size_t index) const
    /// \brief Returns the contents of the \p index cell in the row
    /// \param[in] index The cell number to read
    /// \return The contents in the \p index cell number
    std::string const& operator[](std::size_t index) const
    {
      return m_data[index];
    }
    
    /// \fn std::size_t size() const
    /// \brief Returns the number of cells in #m_data
    /// \return The number of contents
    std::size_t size() const
    {
      return m_data.size();
    }
    
    /// \fn void readNextRow(std::istream & str)
    /// \brief Read and parse the next line from the file \p str
    /// \param[in] str The file to read the next line from
    void readNextRow(std::istream& str)
    {
      /// Get the next line
      std::string         line;
      std::getline(str, line);

      /// Store the line into a std::stringstream
      std::stringstream   lineStream(line);
      std::string         cell;

      /// loop over the stringstream stopping at each ","
      /// add the contents to the #m_data container
      /// first clear the container
      m_data.clear();
      while(std::getline(lineStream, cell, ','))
      {
        m_data.push_back(cell);
      }
      // This checks for a trailing comma with no data after it.
      if (!lineStream && cell.empty())
      {
        // If there was a trailing comma then add an empty element.
        m_data.push_back("");
      }
    } // end void readNextRow

  private:
    std::vector<std::string>    m_data; ///< contents of the csv file for a given line
};

/*!**********************************************
 * \fn std::istream & operator>>(std::istream & str, CSVRow & data)
 * \brief Override the operator>> read a file and save it to a CSV object
 * \param[in] str The input file
 * \param[in,out] data Object to the #CSVRow class to store the next line to
 * \return Return the pointer to \p str
 ***********************************************/
std::istream& operator>>(std::istream& str, CSVRow& data)
{
  data.readNextRow(str);
  return str;
}

/*!**********************************************
 * \class CSVIterator
 * \brief Class to handle looping over the file in an iterator style
 ***********************************************/
class CSVIterator
{   
  public:
    typedef std::input_iterator_tag     iterator_category;
    typedef CSVRow                      value_type;
    typedef std::size_t                 difference_type;
    typedef CSVRow*                     pointer;
    typedef CSVRow&                     reference;

    CSVIterator(std::istream& str)  :m_str(str.good()?&str:NULL) { ++(*this); }
    CSVIterator()                   :m_str(NULL) {}

    // Pre Increment
    CSVIterator& operator++()               {if (m_str) { if (!((*m_str) >> m_row)){m_str = NULL;}}return *this;}
    // Post increment
    CSVIterator operator++(int)             {CSVIterator    tmp(*this);++(*this);return tmp;}
    CSVRow const& operator*()   const       {return m_row;}
    CSVRow const* operator->()  const       {return &m_row;}

    bool operator==(CSVIterator const& rhs) {return ((this == &rhs) || ((this->m_str == NULL) && (rhs.m_str == NULL)));}
    bool operator!=(CSVIterator const& rhs) {return !((*this) == rhs);}
  private:
    std::istream*       m_str;
    CSVRow              m_row;
};


#endif // #ifndef CSVrow_h
