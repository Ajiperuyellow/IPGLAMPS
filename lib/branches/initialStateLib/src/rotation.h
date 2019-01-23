//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/rotation.h $
//$LastChangedDate: 2015-03-13 08:13:52 +0100 (Fr, 13. Mär 2015) $
//$LastChangedRevision: 2123 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief Declaration for routines related to rotation of 4 vectors
 * 
 * This is very closely related to the implementation of class
 * lorentz. Some code may be duplicated. Here a common implementation
 * of something like matrix4D as in vector4D would be very much
 * appreciated. Then boosting and rotating is only applying the matrx
 * multiplications.  
 */

#ifndef ROTATION_H
#define ROTATION_H

#include "bampsvector.h"

#include <iostream>
#include <iomanip>
#include <stdexcept>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


/** 
 * @brief class to implement rotations for vector4D
 * 
 * The internal representation is independent of the user interface.
 **/
class rotation
{
public:
  /** 
   * @brief The type of the numerical representation
   **/
  typedef double Scalar;

protected:
  /** 
   * @brief The internal array to store the values
   **/
  SSE_ALIGNED(Scalar) Arr[16];
  
public:
  
  /** 
   * @brief Constructor
   * 
   * Default constructor of a zero rotation
   */
  rotation( ) 
  { 
    SSE_ALIGNED(static const Scalar) Arr1[16] = 
      { 1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0 };
    for (unsigned int i=0;i<16;i++) Arr[i] = Arr1[i];
  };

  ///// ROTATE ONE VECTOR /////

  /** 
   * @brief Return a rotated vector
   **/
  vector4D rotate(const vector4D & x) const;

  /** 
   * @brief Return a rotated-back vector
   **/
  vector4D rotateInv(const vector4D & x) const;

  /** 
   * @brief Operator overload: Return a rotated vector
   *
   * The multiplication of a rotation with a vector returns a rotated vector.
   **/
  vector4D operator* (const vector4D & x) const 
  { return rotate(x); }

  /** 
   * @brief Rotate one vector
   **/
  void rotate(const vector4D & x, vector4D & xNew) const 
  { xNew = rotate(x); }; 

  /** 
   * @brief Rotate-back one vector
   **/
  void rotateInv(const vector4D & x, vector4D & xNew) const 
  { xNew = rotateInv(x); }; 

  /** 
   * @brief Output routine
   *
   * Writes out the transformation matrix in human readable form.
   **/
  friend std::ostream& operator<<(std::ostream &os, const rotation &obj)
  {
    os << std::setprecision(15) << " (" << obj.Arr[0] 
       << ","  << obj.Arr[3] 
       << ","  << obj.Arr[2] 
       << ","  << obj.Arr[1] << ") " << std::endl;
    for (unsigned int i=3;i>0;i--)
      os << std::setprecision(15) << " (" << obj.Arr[0+4*i] 
         << ","  << obj.Arr[3+4*i] 
         << ","  << obj.Arr[2+4*i] 
         << ","  << obj.Arr[1+4*i] << ") " << std::endl;
    
    return os;
  }
  
  FREE_STORE_OPERATORS_ALIGNED
} __attribute__((aligned(VectorAlignment)));


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


/** 
 * @brief class to implement rotations for vector4D around z-axis
 **/
class rotationZ: public rotation
{
public:
  
  /** 
   * @brief Constructor
   * 
   * Default constructor of a zero rotation
   */
  rotationZ( ) :
    rotation()
  { };

  /** 
   * @brief Constructor with given angle
   **/
  rotationZ( const double angle ) :
    rotation()
  {
    const double c = cos(angle);
    const double s = sin(angle);

    // X: Mem[3], Y: Mem[2]

    Arr[2+4*2] = Arr[3+4*3] = c;
    Arr[2+4*3] = -s;
    Arr[3+4*2] =  s;
  }


};


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif
