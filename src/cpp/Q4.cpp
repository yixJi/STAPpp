/******************************************************************************
 *  STAP++ – CQ4 element implementation (plane stress)                        *
 *  Key points                                                                *
 *  1. Material matrix D does NOT include thickness t                         *
 *  2. Thickness t is applied once per Gauss point                            *
 *  3. Stresses are evaluated at element centre (ξ = η = 0)                   *
 ******************************************************************************/

#include "Q4.h"
#include <cstring>          // std::memset
#include <cstdlib>          // std::exit
#include <iomanip>
#include <iostream>

using namespace std;

/* ---- constructor --------------------------------------------------------- */
CQ4::CQ4()
{
    NEN_  = 4;                      /* nodes per element */
    nodes_ = new CNode*[4];

    ND_   = 8;                      /* element DOF = 4 × 2 */
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

/* ---- element input ------------------------------------------------------- */
bool CQ4::Read(ifstream& input,
               CMaterial* materialSets,
               CNode*     nodeList)
{
    unsigned n1, n2, n3, n4, mset;
    input >> n1 >> n2 >> n3 >> n4 >> mset;

    nodes_[0] = &nodeList[n1 - 1];
    nodes_[1] = &nodeList[n2 - 1];
    nodes_[2] = &nodeList[n3 - 1];
    nodes_[3] = &nodeList[n4 - 1];

    ElementMaterial_ = materialSets + mset - 1;
    return true;
}

/* ---- element output (for .out file) ------------------------------------- */
void CQ4::Write(COutputter& out)
{
    out << setw(9) << nodes_[0]->NodeNumber
        << setw(9) << nodes_[1]->NodeNumber
        << setw(9) << nodes_[2]->NodeNumber
        << setw(9) << nodes_[3]->NodeNumber
        << setw(9) << ElementMaterial_->nset << endl;
}

/* ---- shape functions & local derivatives -------------------------------- */
void CQ4::ShapeFunc(double xi,  double eta,
                    double N[4],
                    double dNdxi[4],
                    double dNdeta[4]) const
{
    N[0] = 0.25 * (1 - xi) * (1 - eta);
    N[1] = 0.25 * (1 + xi) * (1 - eta);
    N[2] = 0.25 * (1 + xi) * (1 + eta);
    N[3] = 0.25 * (1 - xi) * (1 + eta);

    dNdxi [0] = -0.25 * (1 - eta);
    dNdxi [1] =  0.25 * (1 - eta);
    dNdxi [2] =  0.25 * (1 + eta);
    dNdxi [3] = -0.25 * (1 + eta);

    dNdeta[0] = -0.25 * (1 - xi);
    dNdeta[1] = -0.25 * (1 + xi);
    dNdeta[2] =  0.25 * (1 + xi);
    dNdeta[3] =  0.25 * (1 - xi);
}

/* ---- element stiffness (lower-tri. packed) ------------------------------ */
void CQ4::ElementStiffness(double* Matrix)
{
    std::memset(Matrix, 0, sizeof(double) * ND_ * (ND_ + 1) / 2);

    auto* mat = static_cast<CQ4Material*>(ElementMaterial_);
    const double E  = mat->E;
    const double nu = mat->nu;
    const double t  = mat->thickness;

    const double coef = E / (1.0 - nu * nu);
    const double D[3][3] = {
        { 1.0      ,  nu       , 0.0 },
        { nu       ,  1.0      , 0.0 },
        { 0.0      ,  0.0      , (1.0 - nu) / 2.0 }
    };

    /* 2×2 Gauss integration */
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
        {
            const double xi  = gp [a];
            const double eta = gp [b];

            /* N and local derivatives */
            double N[4], dNdxi[4], dNdeta[4];
            ShapeFunc(xi, eta, N, dNdxi, dNdeta);

            /* Jacobian */
            double J[2][2] = {{0}};
            for (int k = 0; k < 4; ++k)
            {
                J[0][0] += dNdxi [k] * nodes_[k]->XYZ[0];
                J[0][1] += dNdxi [k] * nodes_[k]->XYZ[1];
                J[1][0] += dNdeta[k] * nodes_[k]->XYZ[0];
                J[1][1] += dNdeta[k] * nodes_[k]->XYZ[1];
            }
            const double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
            if (detJ <= 0.0)
            {
                cerr << "*** Error *** CQ4 element has non-positive Jacobian.\n";
                std::exit(EXIT_FAILURE);
            }

            /* inverse J */
            const double invJ[2][2] = {
                {  J[1][1]/detJ , -J[0][1]/detJ },
                { -J[1][0]/detJ ,  J[0][0]/detJ }
            };

            /* derivatives wrt x,y */
            double dNdx[4], dNdy[4];
            for (int k = 0; k < 4; ++k)
            {
                dNdx[k] = invJ[0][0]*dNdxi[k] + invJ[0][1]*dNdeta[k];
                dNdy[k] = invJ[1][0]*dNdxi[k] + invJ[1][1]*dNdeta[k];
            }

            /* B-matrix */
            double B[3][8] = {{0}};
            for (int k = 0; k < 4; ++k)
            {
                B[0][2*k    ] = dNdx[k];
                B[1][2*k + 1] = dNdy[k];
                B[2][2*k    ] = dNdy[k];
                B[2][2*k + 1] = dNdx[k];
            }

            /* Ke = ∫ Bt D B */
            double Ke[8][8] = {{0}};
            for (int r = 0; r < 8; ++r)
                for (int c = 0; c <= r; ++c)        /* only lower-tri needed */
                    for (int m = 0; m < 3; ++m)
                        for (int n = 0; n < 3; ++n)
                            Ke[r][c] += B[m][r] * coef * D[m][n] * B[n][c];

            const double factor = t * detJ * wgt[a] * wgt[b];

            /* assemble to lower-triangle packed array */
            unsigned id = 0;
            for (int i = 0; i < 8; ++i)
                for (int j = i; j >=0; --j, ++id)
                    Matrix[id] += Ke[i][j] * factor;
        }
}

/* ---- stresses at element centre (ξ=η=0) -------------------------------- */
void CQ4::ElementStress(double* stress, double* disp)
{
    /* derivatives in local space */
    double N[4], dNdxi[4], dNdeta[4];
    ShapeFunc(0.0, 0.0, N, dNdxi, dNdeta);

    double J[2][2] = {{0}};
    for (int k = 0; k < 4; ++k)
    {
        J[0][0] += dNdxi [k] * nodes_[k]->XYZ[0];
        J[0][1] += dNdxi [k] * nodes_[k]->XYZ[1];
        J[1][0] += dNdeta[k] * nodes_[k]->XYZ[0];
        J[1][1] += dNdeta[k] * nodes_[k]->XYZ[1];
    }
    const double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];

    const double invJ[2][2] = {
        {  J[1][1]/detJ , -J[0][1]/detJ },
        { -J[1][0]/detJ ,  J[0][0]/detJ }
    };

    double dNdx[4], dNdy[4];
    for (int k = 0; k < 4; ++k)
    {
        dNdx[k] = invJ[0][0]*dNdxi[k] + invJ[0][1]*dNdeta[k];
        dNdy[k] = invJ[1][0]*dNdxi[k] + invJ[1][1]*dNdeta[k];
    }

    /* B-matrix */
    double B[3][8] = {{0}};
    for (int k = 0; k < 4; ++k)
    {
        B[0][2*k    ] = dNdx[k];
        B[1][2*k + 1] = dNdy[k];
        B[2][2*k    ] = dNdy[k];
        B[2][2*k + 1] = dNdx[k];
    }

    /* element displacement vector */
    double u[8] = {0};
    for (int i = 0; i < 8; ++i)
        if (LocationMatrix_[i])
            u[i] = disp[LocationMatrix_[i] - 1];

    double strain[3] = {0};
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 8; ++k)
            strain[i] += B[i][k] * u[k];

    auto* mat = static_cast<CQ4Material*>(ElementMaterial_);
    const double E  = mat->E;
    const double nu = mat->nu;
    const double coef = E / (1.0 - nu * nu);

    const double D[3][3] = {
        { 1.0,  nu, 0.0 },
        { nu ,  1.0, 0.0 },
        { 0.0,  0.0, (1.0 - nu) / 2.0 }
    };

    for (int i = 0; i < 3; ++i)
    {
        stress[i] = 0.0;
        for (int j = 0; j < 3; ++j)
            stress[i] += coef * D[i][j] * strain[j];
    }
}

/* ---- location matrix ---------------------------------------------------- */
void CQ4::GenerateLocationMatrix()
{
    unsigned int idx = 0;
    for (unsigned n = 0; n < NEN_; ++n)
        for (unsigned d = 0; d < 2; ++d)
            LocationMatrix_[idx++] = nodes_[n]->bcode[d];
}
