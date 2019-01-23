//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/tools.h $
//$LastChangedDate: 2015-06-10 23:24:51 +0200 (Mi, 10. Jun 2015) $
//$LastChangedRevision: 2176 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file
 * @brief Some useful routines
 **/

#ifndef TOOLS_H
#define TOOLS_H

#include <cstddef>
#include <fstream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>


namespace ns_casc
{
  /** 
   * @brief check whether file is empty
   * 
   * side effect: the file position is set to the end of the file.
   * @param file The file to check
   * 
   * @return true if file length is zero, i.e. file is empty
   */
  inline bool is_empty(std::fstream& file)
  {
    file.seekp(0, std::ios::end);  
    return (file.tellp() == 0);
  }

  /** 
   * @brief check whether file exists
   * 
   * @param name The name of the file to check
   * 
   * @return true if the given file exits (and is readable)
   */
  inline bool exists_file(const std::string& name) 
  {
    if (FILE *file = fopen(name.c_str(), "r")) 
    {
      fclose(file);
      return true;
    } 
    else 
    {
      return false;
    }   
  }


  /**
   * @brief Helper routine to set all entries of a boost::multi_array to
   * zero
   * 
   * @tparam T The type of the array
   * @param array The array to zero
   **/
  template <class T>
  inline void SetZeroArray(T & array)
  {
    for (auto it=array.data(); it!=(array.data()+array.num_elements()); ++it)
    {
      (*it) = 0.0;
    }
  }

  /**
   * @brief Report the (used) size of some std::vector
   */
  template <class T>
  inline std::size_t BytesUsed(std::vector<T> & vec)
  {
    return sizeof(T)*vec.size();
  }

  /**
   * @brief Report the (usable) size of some std::vector
   */
  template <class T>
  inline std::size_t BytesAllocated(std::vector<T> & vec)
  {
    return sizeof(T)*vec.size();
  }

}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @brief Helper class for initializing std::map container
 *
 *
 * @tparam T The type of the first entry of the map
 * @tparam U The type of the second entry of the map
 *
 * Usage: 
 * @code
 * std::map mymap = create_map<int, int >(1,2)(3,4)(5,6);
 * @endcode
 */
template <typename T, typename U>
class create_map
{
private:
  std::map<T, U> m_map;
  
public:
  create_map( const T& key, const U& val )
  {
    m_map[key] = val;
  }

  create_map<T, U>& operator()( const T& key, const U& val )
  {
    m_map[key] = val;
    return *this;
  }

  operator std::map<T, U>()
  {
    return m_map;
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @brief Helper class for initializing std::unordered_map container
 *
 * @tparam T The type of the first entry of the map
 * @tparam U The type of the second entry of the map
 *
 * Usage: 
 * @code
 * std::unordered map mymap = create_unordered_map<int, int>(1,2)(3,4)(5,6);
 * @endcode
 */
template <typename T, typename U>
class create_unordered_map
{
private:
  std::unordered_map<T, U> m_unordered_map;
  
public:
  create_unordered_map( const T& key, const U& val )
  {
    m_unordered_map[key] = val;
  }

  create_unordered_map<T, U>& operator()( const T& key, const U& val )
  {
    m_unordered_map[key] = val;
    return *this;
  }

  operator std::unordered_map<T, U>()
  {
    return m_unordered_map;
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @brief Helper class for initializing std::vector container
 *
 * @tparam T The type of entries of the vector
 *
 * Usage: 
 * @code
 * std::vector myvec = create_vector<int>(1)(2)(3);
 * @endcode
 */
template <typename T>
class create_vector
{
private:
  std::vector<T> m_vec;
  
public:
  create_vector( const T& val )
  {
    m_vec.push_back( val );
  }

  create_vector<T>& operator()( const T& val )
  {
    m_vec.push_back( val );
    return *this;
  }

  operator std::vector<T>()
  {
    return m_vec;
  }
};


#endif
