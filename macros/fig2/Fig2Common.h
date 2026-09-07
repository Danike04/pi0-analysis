/*
 * Fig2Common.h
 *
 * Shared kinematic definitions and constants used by the later Fig. 2
 * coherent/incoherent analysis macros. In particular it defines the four
 * photon-energy intervals, the 8/10/12 degree DeltaPhi cuts, pi0 and 4He
 * masses, photon four-vectors, and the CM missing-energy/opening-angle
 * kinematics used by the template fillers.
 *
 * Keep this header in the same directory (or include path) as the macros that
 * include it. Analysis definitions below are preserved from the original.
 */
#ifndef FIG2_COMMON_H
#define FIG2_COMMON_H

#include <cmath>
#include <algorithm>
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

namespace Fig2 {

static const double kPi0MassMeV = 134.9768;
static const double kHe4MassMeV = 3727.3794;

static const int kNEnergyBins = 4;
static const double kEnergyCenter[kNEnergyBins] = {224.0, 294.0, 320.0, 366.0};

/*
 * Exact photon-energy intervals used in the existing data/coherent plots.
 * The displayed labels are the historical representative energies used in
 * Fig. 2, while these are the actual accepted intervals.
 */
static const double kEnergyLow[kNEnergyBins]  = {223.0, 283.0, 319.0, 356.0};
static const double kEnergyHigh[kNEnergyBins] = {234.0, 294.0, 330.0, 366.0};

static const int kNCuts = 3;
static const double kDeltaPhiCutDeg[kNCuts] = {8.0, 10.0, 12.0};

inline int FindEnergyBin(double eGammaMeV)
{
    for (int i = 0; i < kNEnergyBins; ++i) {
        if (eGammaMeV >= kEnergyLow[i] && eGammaMeV <= kEnergyHigh[i])
            return i;
    }
    return -1;
}

inline double DegToRad(double x) { return x * TMath::DegToRad(); }
inline double RadToDeg(double x) { return x * TMath::RadToDeg(); }

inline double WrapPhi(double phi)
{
    while (phi >  TMath::Pi()) phi -= 2.0*TMath::Pi();
    while (phi <= -TMath::Pi()) phi += 2.0*TMath::Pi();
    return phi;
}

inline TVector3 PhotonDirection(double thetaRad, double phiRad)
{
    return TVector3(std::sin(thetaRad)*std::cos(phiRad),
                    std::sin(thetaRad)*std::sin(phiRad),
                    std::cos(thetaRad));
}

inline TLorentzVector PhotonP4(double energyMeV, double thetaRad, double phiRad)
{
    const TVector3 u = PhotonDirection(thetaRad, phiRad);
    return TLorentzVector(energyMeV*u.X(),
                          energyMeV*u.Y(),
                          energyMeV*u.Z(),
                          energyMeV);
}

/* Paper Eq. (1). */
inline double InvariantMass(double e1, double e2, double openingAngleRad)
{
    const double m2 = 2.0*e1*e2*(1.0 - std::cos(openingAngleRad));
    return m2 > 0.0 ? std::sqrt(m2) : 0.0;
}

/* Paper Eq. (3): expected pion CM energy under coherent kinematics. */
inline double ExpectedPionCmEnergy(double eGammaMeV,
                                   double mPiMeV = kPi0MassMeV,
                                   double mTargetMeV = kHe4MassMeV)
{
    const double s = mTargetMeV*mTargetMeV + 2.0*eGammaMeV*mTargetMeV;
    if (s <= 0.0) return 0.0;
    return (s + mPiMeV*mPiMeV - mTargetMeV*mTargetMeV)/(2.0*std::sqrt(s));
}

inline double CmBeta(double eGammaMeV, double mTargetMeV = kHe4MassMeV)
{
    return eGammaMeV/(eGammaMeV + mTargetMeV);
}

inline double LorentzGammaFromBeta(double beta)
{
    const double oneMinusBeta2 = 1.0 - beta*beta;
    return oneMinusBeta2 > 0.0 ? 1.0/std::sqrt(oneMinusBeta2) : 0.0;
}

/* Paper Eq. (4), with the standard Lorentz gamma = 1/sqrt(1-beta^2). */
inline double MeasuredPionCmEnergy(double eGammaMeV,
                                  const TLorentzVector& pi0Lab,
                                  double mTargetMeV = kHe4MassMeV)
{
    const double beta = CmBeta(eGammaMeV, mTargetMeV);
    const double gamma = LorentzGammaFromBeta(beta);
    return gamma*(pi0Lab.E() - beta*pi0Lab.Pz());
}

/* Paper Eq. (5). Negative values form the incoherent/break-up side. */
inline double DeltaECoherentHypothesis(double eGammaMeV,
                                      const TLorentzVector& pi0Lab,
                                      double mPiMeV = kPi0MassMeV,
                                      double mTargetMeV = kHe4MassMeV)
{
    return MeasuredPionCmEnergy(eGammaMeV, pi0Lab, mTargetMeV)
         - ExpectedPionCmEnergy(eGammaMeV, mPiMeV, mTargetMeV);
}

/*
 * Solve gamma + A -> pi0 + A for the pion LAB energy at a measured pion
 * polar angle. This supplies the coherent reference energy used in Phi_min.
 */
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

    double disc = A*A*C*C - D*(C*C + B*B*mPiMeV*mPiMeV);
    if (disc < 0.0 && disc > -1e-7) disc = 0.0;
    if (disc < 0.0) return -1.0;

    const double root = std::sqrt(disc);
    const double candidates[2] = {(A*C + root)/D, (A*C - root)/D};

    double bestE = -1.0;
    double bestResidual = 1e99;

    for (int i = 0; i < 2; ++i) {
        const double E = candidates[i];
        if (E < mPiMeV || E > eGammaMeV + mPiMeV + 50.0) continue;
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

/* Paper Eq. (2). */
inline double MinimumOpeningAngle(double pionLabEnergyMeV,
                                  double mPiMeV = kPi0MassMeV)
{
    if (pionLabEnergyMeV <= mPiMeV) return TMath::Pi();
    double ratio = std::sqrt(std::max(0.0,
                      pionLabEnergyMeV*pionLabEnergyMeV - mPiMeV*mPiMeV))
                 / pionLabEnergyMeV;
    ratio = std::max(-1.0, std::min(1.0, ratio));
    return 2.0*std::acos(ratio);
}

/*
 * DeltaPhi = measured gamma-gamma opening angle minus the minimum opening
 * angle expected for a coherently produced pi0 at the reconstructed pi0 angle.
 */
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
        CoherentPionLabEnergy(eGammaMeV, pi0.Theta(), mPiMeV, mTargetMeV);
    if (expectedLabE <= mPiMeV) return 1e9;

    const double phiMin = MinimumOpeningAngle(expectedLabE, mPiMeV);
    return RadToDeg(opening - phiMin);
}

} // namespace Fig2

#endif
