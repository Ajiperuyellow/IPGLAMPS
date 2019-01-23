//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/pdfinterface.h $
//$LastChangedDate: 2012-03-21 22:26:37 +0100 (水, 21  3月 2012) $
//$LastChangedRevision: 520 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef DOUBLEPDF_H
#define DOUBLEPDF_H

#include <stdexcept>
#include <string>



/**
 * @brief Class to calculate the pdf's for beam and target nucleon
 *
 * Here we calculate the pdf's for the beam and target nucleon
 * simultanously.
 */
class interfacePDF_generic
{
public:
  /** 
   * @brief Constructor
   */
    interfacePDF_generic( ) {};

  /** 
   * @brief The actual routine
   *
   * @param[in] Q2 The scale (in GeV^2).
   * @param[in] x1,x2 The x values.
   * @param[in] nuc1,nuc2 The charges of the target nucleons: 
   *   =1 for proton, =0 for neutron
   * @param[out] F1,F2 The pdf values. 
   *
   * The indices correspond to:
   *   0->g, 1->u,2->ubar, 3->d,4->dbar, 5->s,6->sbar, 7->c,8->cbar,
   *   9->b,10->bbar, 11->t,12->tbar
   */
  virtual void eval(double Q2, double x1, double x2, int nuc1, int nuc2, double F1[13], double F2[13]) = 0;

};



/**
 * @brief Class to provide GRV pdf's
 *
 * This class uses the 'old' Gluck-Reya-Vogt parametrizations implemented by Xu and
 * should be superseded by the corresponding LHAPDF implementation.
 *
 * Reference: Z.Phys C 67, 433-447(1995)
 */
class interfacePDF_GRV : public interfacePDF_generic
{
public:
  interfacePDF_GRV( ) {};

  virtual void eval(double Q2, double x1, double x2, int nuc1, int nuc2, double F1[13], double F2[13]);
};



/**
 * @brief Class to provide pdf's implemented in LHAPDF
 *
 * This class is a wrapper around the corresponding LHAPDF routines 
 *
 * By giving a valid name to the constructor, the corresponding PDF is
 * used. Usage of nPDFs, e.g. EPS09, is also possible.
 *
 * The LHAPDF routines are only available, if the
 * corresponding library is found by the build system. Please note,
 * that even if the library is found, some data files may be missing
 * and have to be installed by you by hand.
 **/
class interfacePDF_LHAPDF : public interfacePDF_generic
{
public:
  interfacePDF_LHAPDF ( const std::string& name, const unsigned int member, const bool useGrid, const bool _nPDF = false, const std::string _nPDFname = "", int _A1 = 1, int _A2 = 1 );

  virtual void eval(double Q2, double x1, double x2, int nuc1, int nuc2, double F1[13], double F2[13]);

private:
  int A1, A2;
  bool useNuclearPDFs;
  
protected:
  /**
   * @brief Map the LHAPDF values to the internal encoding.
   *
   * LHAPDF:
   *  - 0: tbar
   *  - 1: bbar
   *  - 2: cbar
   *  - 3: sbar
   *  - 4: ubar
   *  - 5: dbar
   *  - 6: gluon
   *  - 7: d
   *  - 8: u
   *  - 9: s
   *  - 10: c
   *  - 11: b
   *  - 12: t
   *
   * Own encoding:
   *  -     0: gluon
   *  -  1, 2: u, ubar
   *  -  3, 4: d, dbar
   *  -  5, 6: s, sbar
   *  -  7, 8: c, cbar
   *  -  9,10: b, bbar
   *  - 11,12: t, tbar
   *
   * Since the LHAPDF are for proton only (charge==1) , we also
   * exchange the u,ubar contributions with the d,dbar contributions,
   * if the the neutron pdf is to be calculated (charge==0).
   */
  void Map(const int charge, const double F_LHAPDF[13], double F[13]);
};



/** 
 * @brief exception class for handling unexpected critical behaviour
 * within generation of pdf
 **/
class ePDF_error : public std::runtime_error
{
public:
  explicit ePDF_error(const std::string& what) : std::runtime_error(what) {};
  virtual ~ePDF_error() throw() {};
};



#endif