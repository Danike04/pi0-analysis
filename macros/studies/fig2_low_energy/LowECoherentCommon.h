/*
 * LowECoherentCommon.h
 *
 * Shared constants and kinematic helper functions for the dedicated low-energy
 * coherent Fig. 2 studies.  It defines the 135-215 MeV energy bins, the
 * 8/10/12-degree DeltaPhi cuts, coherent two-body kinematics, CM missing energy
 * and the minimum gamma-gamma opening-angle quantity used by that branch.
 *
 * This header is preserved functionally unchanged from the supplied analysis.
 */
#ifndef LOWE_COHERENT_COMMON_H
#define LOWE_COHERENT_COMMON_H

#include <cmath>
#include <algorithm>
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

namespace LowECoh {

static const double kPi0MassMeV = 134.9768;
static const double kHe4MassMeV = 3727.3794;

static const int kNEnergyBins = 4;
static const double kEnergyLow[kNEnergyBins]  = {135.0, 155.0, 175.0, 195.0};
static const double kEnergyHigh[kNEnergyBins] = {155.0, 175.0, 195.0, 215.0};
static const double kEnergyCenter[kNEnergyBins] = {145.0, 165.0, 185.0, 205.0};

static const int kNCuts = 3;
static const double kDeltaPhiCutDeg[kNCuts] = {8.0, 10.0, 12.0};

inline int FindEnergyBin(double eGammaMeV)
{
    for (int i = 0; i < kNEnergyBins; ++i) {
        const bool inBin =
            eGammaMeV >= kEnergyLow[i] &&
            (eGammaMeV < kEnergyHigh[i] ||
             (i == kNEnergyBins - 1 && eGammaMeV <= kEnergyHigh[i]));
        if (inBin) return i;
    }
    return -1;
}

inline double DegToRad(double x) { return x*TMath::DegToRad(); }
inline double RadToDeg(double x) { return x*TMath::RadToDeg(); }

inline TVector3 PhotonDirection(double thetaRad, double phiRad)
{
    return TVector3(std::sin(thetaRad)*std::cos(phiRad),
                    std::sin(thetaRad)*std::sin(phiRad),
                    std::cos(thetaRad));
}

inline TLorentzVector PhotonP4(double energyMeV,
                              double thetaRad,
                              double phiRad)
{
    const TVector3 u = PhotonDirection(thetaRad, phiRad);
    return TLorentzVector(energyMeV*u.X(),
                          energyMeV*u.Y(),
                          energyMeV*u.Z(),
                          energyMeV);
}

inline double ExpectedPionCmEnergy(double eGammaMeV,
                                   double mPiMeV = kPi0MassMeV,
                                   double mTargetMeV = kHe4MassMeV)
{
    const double s = mTargetMeV*mTargetMeV +
                     2.0*eGammaMeV*mTargetMeV;
    if (s <= 0.0) return 0.0;
    return (s + mPiMeV*mPiMeV - mTargetMeV*mTargetMeV)/
           (2.0*std::sqrt(s));
}

inline double CmBeta(double eGammaMeV,
                     double mTargetMeV = kHe4MassMeV)
{
    return eGammaMeV/(eGammaMeV + mTargetMeV);
}

inline double LorentzGammaFromBeta(double beta)
{
    const double x = 1.0 - beta*beta;
    return x > 0.0 ? 1.0/std::sqrt(x) : 0.0;
}

inline double MeasuredPionCmEnergy(double eGammaMeV,
                                  const TLorentzVector& pi0Lab,
                                  double mTargetMeV = kHe4MassMeV)
{
    const double beta = CmBeta(eGammaMeV, mTargetMeV);
    const double gamma = LorentzGammaFromBeta(beta);
    return gamma*(pi0Lab.E() - beta*pi0Lab.Pz());
}

inline double DeltaECoherentHypothesis(double eGammaMeV,
                                      const TLorentzVector& pi0Lab,
                                      double mPiMeV = kPi0MassMeV,
                                      double mTargetMeV = kHe4MassMeV)
{
    return MeasuredPionCmEnergy(eGammaMeV, pi0Lab, mTargetMeV)
         - ExpectedPionCmEnergy(eGammaMeV, mPiMeV, mTargetMeV);
}

inline double CoherentPionLabEnergy(double eGammaMeV,
                                    double thetaPiLabRad,
                                    double mPiMeV = kPi0MassMeV,
                                    double mTargetMeV = kHe4MassMeV)
{
    const double k = eGammaMeV;
    const double c = std::cos(thetaPiLabRad);
    const double A = k + mTargetMeV;
    const double B = k*c;
    const double C = k*mTargetMeV + 0.5*mPiMeV*mPiMeV;
    const double D = A*A - B*B;

    if (D <= 0.0) return -1.0;

    double disc = A*A*C*C -
                  D*(C*C + B*B*mPiMeV*mPiMeV);
    if (disc < 0.0 && disc > -1e-7) disc = 0.0;
    if (disc < 0.0) return -1.0;

    const double root = std::sqrt(disc);
    const double candidates[2] =
        {(A*C + root)/D, (A*C - root)/D};

    double bestE = -1.0;
    double bestResidual = 1e99;

    for (int i = 0; i < 2; ++i) {
        const double E = candidates[i];
        if (E < mPiMeV || E > eGammaMeV + mPiMeV + 50.0)
            continue;

        const double p2 = E*E - mPiMeV*mPiMeV;
        if (p2 < 0.0) continue;

        const double p = std::sqrt(p2);
        const double residual = std::fabs(A*E - B*p - C);
        if (residual < bestResidual) {
            bestResidual = residual;
            bestE = E;
        }
    }
    return bestE;
}

inline double MinimumOpeningAngle(double pionLabEnergyMeV,
                                  double mPiMeV = kPi0MassMeV)
{
    if (pionLabEnergyMeV <= mPiMeV)
        return TMath::Pi();

    double ratio =
        std::sqrt(std::max(0.0,
                  pionLabEnergyMeV*pionLabEnergyMeV -
                  mPiMeV*mPiMeV))/pionLabEnergyMeV;

    ratio = std::max(-1.0, std::min(1.0, ratio));
    return 2.0*std::acos(ratio);
}

inline double DeltaPhiCoherentDeg(double eGammaMeV,
                                  const TLorentzVector& gamma1,
                                  const TLorentzVector& gamma2,
                                  double mPiMeV = kPi0MassMeV,
                                  double mTargetMeV = kHe4MassMeV)
{
    const TLorentzVector pi0 = gamma1 + gamma2;
    if (pi0.P() <= 0.0) return 1e9;

    const double opening = gamma1.Vect().Angle(gamma2.Vect());
    const double expectedLabE =
        CoherentPionLabEnergy(eGammaMeV, pi0.Theta(),
                              mPiMeV, mTargetMeV);

    if (expectedLabE <= mPiMeV) return 1e9;

    const double phiMin =
        MinimumOpeningAngle(expectedLabE, mPiMeV);

    return RadToDeg(opening - phiMin);
}

} // namespace LowECoh

#endif
