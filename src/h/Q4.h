#pragma once
/***********************************************************************
 *  CQ4 : 4-node isoparametric plane-stress element (2 DOF / node)     *
 *  ------------------------------------------------------------------ *
 *  Nodal DOF order in every element: { u1 v1 u2 v2 u3 v3 u4 v4 }      *
 ***********************************************************************/

#include "Element.h"
#include "Material.h"

class CQ4 : public CElement
{
public:
    CQ4();
    ~CQ4() override = default;

    /* ---- I/O -------------------------------------------------------- */
    bool Read  (std::ifstream& input,
                CMaterial*     materialSets,
                CNode*         nodeList) override;
    void Write (COutputter& output) override;

    /* ---- FEM routines ---------------------------------------------- */
    void ElementStiffness(double* Matrix) override;          /* lower-tri. packed */
    void ElementStress   (double* stress,                    /* {σx,σy,τxy} */
                          double* disp ) override;           /* global displacement */

    void GenerateLocationMatrix() override;

private:
    /* 2-point Gauss abscissae / weights (±1/√3 , 1) */
    static constexpr double gp [2] = { -0.5773502691896258,
                                        0.5773502691896258 };
    static constexpr double wgt[2] = { 1.0, 1.0 };

    /* shape functions and derivatives at (ξ,η) */
    void ShapeFunc(double  xi,   double  eta,
                   double  N[4],
                   double  dNdxi [4],
                   double  dNdeta[4]) const;
};
