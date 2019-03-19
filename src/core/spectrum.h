
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_SPECTRUM_H
#define PBRT_CORE_SPECTRUM_H

// core/spectrum.h*
#include "pbrt.h"
#include "stringprint.h"

namespace pbrt {

// Spectrum Utility Declarations
static const int sampledLambdaStart = 400;
static const int sampledLambdaEnd = 700;
static const int nSpectralSamples = 60;
extern bool SpectrumSamplesSorted(const Float *lambda, const Float *vals,
                                  int n);
extern void SortSpectrumSamples(Float *lambda, Float *vals, int n);
extern Float AverageSpectrumSamples(const Float *lambda, const Float *vals,
                                    int n, Float lambdaStart, Float lambdaEnd);
inline void XYZToRGB(const Float xyz[3], Float rgb[3]) {
    rgb[0] = 3.240479f * xyz[0] - 1.537150f * xyz[1] - 0.498535f * xyz[2];
    rgb[1] = -0.969256f * xyz[0] + 1.875991f * xyz[1] + 0.041556f * xyz[2];
    rgb[2] = 0.055648f * xyz[0] - 0.204043f * xyz[1] + 1.057311f * xyz[2];
}

inline void RGBToXYZ(const Float rgb[3], Float xyz[3]) {
    xyz[0] = 0.412453f * rgb[0] + 0.357580f * rgb[1] + 0.180423f * rgb[2];
    xyz[1] = 0.212671f * rgb[0] + 0.715160f * rgb[1] + 0.072169f * rgb[2];
    xyz[2] = 0.019334f * rgb[0] + 0.119193f * rgb[1] + 0.950227f * rgb[2];
}

enum class SpectrumType { Reflectance, Illuminant };
extern Float InterpolateSpectrumSamples(const Float *lambda, const Float *vals,
                                        int n, Float l);
extern void Blackbody(const Float *lambda, int n, Float T, Float *Le);
extern void BlackbodyNormalized(const Float *lambda, int n, Float T,
                                Float *vals);

// Spectral Data Declarations
static const int nCIESamples = 471;
extern const Float CIE_X[nCIESamples];
extern const Float CIE_Y[nCIESamples];
extern const Float CIE_Z[nCIESamples];
extern const Float CIE_lambda[nCIESamples];
static const Float CIE_Y_integral = 106.856895;
static const int nRGB2SpectSamples = 32;
extern const Float RGB2SpectLambda[nRGB2SpectSamples];
extern const Float RGBRefl2SpectWhite[nRGB2SpectSamples];
extern const Float RGBRefl2SpectCyan[nRGB2SpectSamples];
extern const Float RGBRefl2SpectMagenta[nRGB2SpectSamples];
extern const Float RGBRefl2SpectYellow[nRGB2SpectSamples];
extern const Float RGBRefl2SpectRed[nRGB2SpectSamples];
extern const Float RGBRefl2SpectGreen[nRGB2SpectSamples];
extern const Float RGBRefl2SpectBlue[nRGB2SpectSamples];
extern const Float RGBIllum2SpectWhite[nRGB2SpectSamples];
extern const Float RGBIllum2SpectCyan[nRGB2SpectSamples];
extern const Float RGBIllum2SpectMagenta[nRGB2SpectSamples];
extern const Float RGBIllum2SpectYellow[nRGB2SpectSamples];
extern const Float RGBIllum2SpectRed[nRGB2SpectSamples];
extern const Float RGBIllum2SpectGreen[nRGB2SpectSamples];
extern const Float RGBIllum2SpectBlue[nRGB2SpectSamples];

// Spectrum Declarations
template <int nSpectrumSamples>
class CoefficientSpectrum {
  public:
    // CoefficientSpectrum Public Methods
    CoefficientSpectrum(Float v = 0.f) {
        for (int i = 0; i < nSpectrumSamples; ++i) c[i] = v;
        DCHECK(!HasNaNs());
    }
#if defined(DEBUG) || defined(PBRT_RUNTIME_RESOLUTION_SPECTRUM)
    CoefficientSpectrum(const CoefficientSpectrum &s) {
        DCHECK(!s.HasNaNs());
        for (int i = 0; i < nSpectrumSamples; ++i) c[i] = s.c[i];
    }

    CoefficientSpectrum &operator=(const CoefficientSpectrum &s) {
        DCHECK(!s.HasNaNs());
        for (int i = 0; i < nSpectrumSamples; ++i) c[i] = s.c[i];
        return *this;
    }
#endif  // DEBUG
    void Print(FILE *f) const {
        fprintf(f, "[ ");
        for (int i = 0; i < nSpectrumSamples; ++i) {
            fprintf(f, "%f", c[i]);
            if (i != nSpectrumSamples - 1) fprintf(f, ", ");
        }
        fprintf(f, "]");
    }
    CoefficientSpectrum &operator+=(const CoefficientSpectrum &s2) {
        DCHECK(!s2.HasNaNs());
        for (int i = 0; i < nSpectrumSamples; ++i) c[i] += s2.c[i];
        return *this;
    }
    CoefficientSpectrum operator+(const CoefficientSpectrum &s2) const {
        DCHECK(!s2.HasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] += s2.c[i];
        return ret;
    }
    CoefficientSpectrum operator-(const CoefficientSpectrum &s2) const {
        DCHECK(!s2.HasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] -= s2.c[i];
        return ret;
    }
    CoefficientSpectrum operator/(const CoefficientSpectrum &s2) const {
        DCHECK(!s2.HasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) {
          CHECK_NE(s2.c[i], 0);
          ret.c[i] /= s2.c[i];
        }
        return ret;
    }
    CoefficientSpectrum operator*(const CoefficientSpectrum &sp) const {
        DCHECK(!sp.HasNaNs());
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] *= sp.c[i];
        return ret;
    }
    CoefficientSpectrum &operator*=(const CoefficientSpectrum &sp) {
        DCHECK(!sp.HasNaNs());
        for (int i = 0; i < nSpectrumSamples; ++i) c[i] *= sp.c[i];
        return *this;
    }
    CoefficientSpectrum operator*(Float a) const {
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] *= a;
        DCHECK(!ret.HasNaNs());
        return ret;
    }
    CoefficientSpectrum &operator*=(Float a) {
        for (int i = 0; i < nSpectrumSamples; ++i) c[i] *= a;
        DCHECK(!HasNaNs());
        return *this;
    }
    friend inline CoefficientSpectrum operator*(Float a,
                                                const CoefficientSpectrum &s) {
        DCHECK(!std::isnan(a) && !s.HasNaNs());
        return s * a;
    }
    CoefficientSpectrum operator/(Float a) const {
        CHECK_NE(a, 0);
        DCHECK(!std::isnan(a));
        CoefficientSpectrum ret = *this;
        for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] /= a;
        DCHECK(!ret.HasNaNs());
        return ret;
    }
    CoefficientSpectrum &operator/=(Float a) {
        CHECK_NE(a, 0);
        DCHECK(!std::isnan(a));
        for (int i = 0; i < nSpectrumSamples; ++i) c[i] /= a;
        return *this;
    }
    bool operator==(const CoefficientSpectrum &sp) const {
        for (int i = 0; i < nSpectrumSamples; ++i)
            if (c[i] != sp.c[i]) return false;
        return true;
    }
    bool operator!=(const CoefficientSpectrum &sp) const {
        return !(*this == sp);
    }
    bool IsBlack() const {
        for (int i = 0; i < nSpectrumSamples; ++i)
            if (c[i] != 0.) return false;
        return true;
    }
    friend CoefficientSpectrum Sqrt(const CoefficientSpectrum &s) {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] = std::sqrt(s.c[i]);
        DCHECK(!ret.HasNaNs());
        return ret;
    }
    template <int n>
    friend inline CoefficientSpectrum<n> Pow(const CoefficientSpectrum<n> &s,
                                             Float e);
    CoefficientSpectrum operator-() const {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] = -c[i];
        return ret;
    }
    friend CoefficientSpectrum Exp(const CoefficientSpectrum &s) {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] = std::exp(s.c[i]);
        DCHECK(!ret.HasNaNs());
        return ret;
    }
    friend std::ostream &operator<<(std::ostream &os,
                                    const CoefficientSpectrum &s) {
        return os << s.ToString();
    }
    std::string ToString() const {
        std::string str = "[ ";
        for (int i = 0; i < nSpectrumSamples; ++i) {
            str += StringPrintf("%f", c[i]);
            if (i + 1 < nSpectrumSamples) str += ", ";
        }
        str += " ]";
        return str;
    }
    CoefficientSpectrum Clamp(Float low = 0, Float high = Infinity) const {
        CoefficientSpectrum ret;
        for (int i = 0; i < nSpectrumSamples; ++i)
            ret.c[i] = pbrt::Clamp(c[i], low, high);
        DCHECK(!ret.HasNaNs());
        return ret;
    }
    Float MaxComponentValue() const {
        Float m = c[0];
        for (int i = 1; i < nSpectrumSamples; ++i)
            m = std::max(m, c[i]);
        return m;
    }
    bool HasNaNs() const {
        for (int i = 0; i < nSpectrumSamples; ++i)
            if (std::isnan(c[i])) return true;
        return false;
    }
    bool Write(FILE *f) const {
        for (int i = 0; i < nSpectrumSamples; ++i)
            if (fprintf(f, "%f ", c[i]) < 0) return false;
        return true;
    }
    bool Read(FILE *f) {
        for (int i = 0; i < nSpectrumSamples; ++i) {
            double v;
            if (fscanf(f, "%lf ", &v) != 1) return false;
            c[i] = v;
        }
        return true;
    }
    Float &operator[](int i) {
        DCHECK(i >= 0 && i < nSpectrumSamples);
        return c[i];
    }
    Float operator[](int i) const {
        DCHECK(i >= 0 && i < nSpectrumSamples);
        return c[i];
    }

    // CoefficientSpectrum Public Data
    static const int nSamples = nSpectrumSamples;

  protected:
    // CoefficientSpectrum Protected Data
    Float c[nSpectrumSamples];
};

class SampledSpectrum : public CoefficientSpectrum<nSpectralSamples> {
  public:
    // SampledSpectrum Public Methods
    SampledSpectrum(Float v = 0.f) : CoefficientSpectrum(v) {}
    SampledSpectrum(const CoefficientSpectrum<nSpectralSamples> &v)
        : CoefficientSpectrum<nSpectralSamples>(v) {}
    static SampledSpectrum FromSampled(const Float *lambda, const Float *v,
                                       int n) {
        // Sort samples if unordered, use sorted for returned spectrum
        if (!SpectrumSamplesSorted(lambda, v, n)) {
            std::vector<Float> slambda(&lambda[0], &lambda[n]);
            std::vector<Float> sv(&v[0], &v[n]);
            SortSpectrumSamples(&slambda[0], &sv[0], n);
            return FromSampled(&slambda[0], &sv[0], n);
        }
        SampledSpectrum r;
        for (int i = 0; i < nSpectralSamples; ++i) {
            // Compute average value of given SPD over $i$th sample's range
            Float lambda0 = Lerp(Float(i) / Float(nSpectralSamples),
                                 sampledLambdaStart, sampledLambdaEnd);
            Float lambda1 = Lerp(Float(i + 1) / Float(nSpectralSamples),
                                 sampledLambdaStart, sampledLambdaEnd);
            r.c[i] = AverageSpectrumSamples(lambda, v, n, lambda0, lambda1);
        }
        return r;
    }
    static void Init() {
        // Compute XYZ matching functions for _SampledSpectrum_
        for (int i = 0; i < nSpectralSamples; ++i) {
            Float wl0 = Lerp(Float(i) / Float(nSpectralSamples),
                             sampledLambdaStart, sampledLambdaEnd);
            Float wl1 = Lerp(Float(i + 1) / Float(nSpectralSamples),
                             sampledLambdaStart, sampledLambdaEnd);
            X.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_X, nCIESamples, wl0,
                                            wl1);
            Y.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Y, nCIESamples, wl0,
                                            wl1);
            Z.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Z, nCIESamples, wl0,
                                            wl1);
        }

        // Compute RGB to spectrum functions for _SampledSpectrum_
        for (int i = 0; i < nSpectralSamples; ++i) {
            Float wl0 = Lerp(Float(i) / Float(nSpectralSamples),
                             sampledLambdaStart, sampledLambdaEnd);
            Float wl1 = Lerp(Float(i + 1) / Float(nSpectralSamples),
                             sampledLambdaStart, sampledLambdaEnd);
            rgbRefl2SpectWhite.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectWhite,
                                       nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectCyan.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectCyan,
                                       nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectMagenta.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectMagenta,
                                       nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectYellow.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectYellow,
                                       nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectRed.c[i] = AverageSpectrumSamples(
                RGB2SpectLambda, RGBRefl2SpectRed, nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectGreen.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectGreen,
                                       nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectBlue.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectBlue,
                                       nRGB2SpectSamples, wl0, wl1);

            rgbIllum2SpectWhite.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectWhite,
                                       nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectCyan.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectCyan,
                                       nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectMagenta.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectMagenta,
                                       nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectYellow.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectYellow,
                                       nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectRed.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectRed,
                                       nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectGreen.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectGreen,
                                       nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectBlue.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectBlue,
                                       nRGB2SpectSamples, wl0, wl1);
        }
    }
    void ToXYZ(Float xyz[3]) const {
        xyz[0] = xyz[1] = xyz[2] = 0.f;
        for (int i = 0; i < nSpectralSamples; ++i) {
            xyz[0] += X.c[i] * c[i];
            xyz[1] += Y.c[i] * c[i];
            xyz[2] += Z.c[i] * c[i];
        }
        Float scale = Float(sampledLambdaEnd - sampledLambdaStart) /
                      Float(CIE_Y_integral * nSpectralSamples);
        xyz[0] *= scale;
        xyz[1] *= scale;
        xyz[2] *= scale;
    }
    Float y() const {
        Float yy = 0.f;
        for (int i = 0; i < nSpectralSamples; ++i) yy += Y.c[i] * c[i];
        return yy * Float(sampledLambdaEnd - sampledLambdaStart) /
               Float(CIE_Y_integral * nSpectralSamples);
    }
    void ToRGB(Float rgb[3]) const {
        Float xyz[3];
        ToXYZ(xyz);
        XYZToRGB(xyz, rgb);
    }
    RGBSpectrum ToRGBSpectrum() const;
    static SampledSpectrum FromRGB(
        const Float rgb[3], SpectrumType type = SpectrumType::Illuminant);
    static SampledSpectrum FromXYZ(
        const Float xyz[3], SpectrumType type = SpectrumType::Reflectance) {
        Float rgb[3];
        XYZToRGB(xyz, rgb);
        return FromRGB(rgb, type);
    }
    SampledSpectrum(const RGBSpectrum &r,
                    SpectrumType type = SpectrumType::Reflectance);

  private:
    // SampledSpectrum Private Data
    static SampledSpectrum X, Y, Z;
    static SampledSpectrum rgbRefl2SpectWhite, rgbRefl2SpectCyan;
    static SampledSpectrum rgbRefl2SpectMagenta, rgbRefl2SpectYellow;
    static SampledSpectrum rgbRefl2SpectRed, rgbRefl2SpectGreen;
    static SampledSpectrum rgbRefl2SpectBlue;
    static SampledSpectrum rgbIllum2SpectWhite, rgbIllum2SpectCyan;
    static SampledSpectrum rgbIllum2SpectMagenta, rgbIllum2SpectYellow;
    static SampledSpectrum rgbIllum2SpectRed, rgbIllum2SpectGreen;
    static SampledSpectrum rgbIllum2SpectBlue;
};

class RGBSpectrum : public CoefficientSpectrum<3> {
    using CoefficientSpectrum<3>::c;

  public:
    // RGBSpectrum Public Methods
    RGBSpectrum(Float v = 0.f) : CoefficientSpectrum<3>(v) {}
    RGBSpectrum(const CoefficientSpectrum<3> &v) : CoefficientSpectrum<3>(v) {}
    RGBSpectrum(const RGBSpectrum &s,
                SpectrumType type = SpectrumType::Reflectance) {
        *this = s;
    }
    static RGBSpectrum FromRGB(const Float rgb[3],
                               SpectrumType type = SpectrumType::Reflectance) {
        RGBSpectrum s;
        s.c[0] = rgb[0];
        s.c[1] = rgb[1];
        s.c[2] = rgb[2];
        DCHECK(!s.HasNaNs());
        return s;
    }
    void ToRGB(Float *rgb) const {
        rgb[0] = c[0];
        rgb[1] = c[1];
        rgb[2] = c[2];
    }
    const RGBSpectrum &ToRGBSpectrum() const { return *this; }
    void ToXYZ(Float xyz[3]) const { RGBToXYZ(c, xyz); }
    static RGBSpectrum FromXYZ(const Float xyz[3],
                               SpectrumType type = SpectrumType::Reflectance) {
        RGBSpectrum r;
        XYZToRGB(xyz, r.c);
        return r;
    }
    Float y() const {
        const Float YWeight[3] = {0.212671f, 0.715160f, 0.072169f};
        return YWeight[0] * c[0] + YWeight[1] * c[1] + YWeight[2] * c[2];
    }
    static RGBSpectrum FromSampled(const Float *lambda, const Float *v, int n) {
        // Sort samples if unordered, use sorted for returned spectrum
        if (!SpectrumSamplesSorted(lambda, v, n)) {
            std::vector<Float> slambda(&lambda[0], &lambda[n]);
            std::vector<Float> sv(&v[0], &v[n]);
            SortSpectrumSamples(&slambda[0], &sv[0], n);
            return FromSampled(&slambda[0], &sv[0], n);
        }
        Float xyz[3] = {0, 0, 0};
        for (int i = 0; i < nCIESamples; ++i) {
            Float val = InterpolateSpectrumSamples(lambda, v, n, CIE_lambda[i]);
            xyz[0] += val * CIE_X[i];
            xyz[1] += val * CIE_Y[i];
            xyz[2] += val * CIE_Z[i];
        }
        Float scale = Float(CIE_lambda[nCIESamples - 1] - CIE_lambda[0]) /
                      Float(CIE_Y_integral * nCIESamples);
        xyz[0] *= scale;
        xyz[1] *= scale;
        xyz[2] *= scale;
        return FromXYZ(xyz);
    }
};

// Spectrum Inline Functions
template <int nSpectrumSamples>
inline CoefficientSpectrum<nSpectrumSamples> Pow(
    const CoefficientSpectrum<nSpectrumSamples> &s, Float e) {
    CoefficientSpectrum<nSpectrumSamples> ret;
    for (int i = 0; i < nSpectrumSamples; ++i) ret.c[i] = std::pow(s.c[i], e);
    DCHECK(!ret.HasNaNs());
    return ret;
}

inline RGBSpectrum Lerp(Float t, const RGBSpectrum &s1, const RGBSpectrum &s2) {
    return (1 - t) * s1 + t * s2;
}

inline SampledSpectrum Lerp(Float t, const SampledSpectrum &s1,
                            const SampledSpectrum &s2) {
    return (1 - t) * s1 + t * s2;
}

void ResampleLinearSpectrum(const Float *lambdaIn, const Float *vIn, int nIn,
                            Float lambdaMin, Float lambdaMax, int nOut,
                            Float *vOut);


#ifdef PBRT_RUNTIME_RESOLUTION_SPECTRUM

#include "error.h"

// Adding realtime spectrum resolution
class RuntimeResolutionSpectrum
{
    enum class RuntimeSpectrumType { 
        SampledSpectrumType, 
        RGBSpectrumType, 
        count 
    };
    static RuntimeSpectrumType spectrumType;

    union RuntimeSpectrum
    {
        SampledSpectrum* pSS;
        RGBSpectrum* pTS;
    };
    RuntimeSpectrum runtimeSpectrum;

public:
    
    static void setSampledSpectrum()
    {
        spectrumType = RuntimeSpectrumType::SampledSpectrumType;
        nSamples = 60;
    }

    static void setRGBSpectrum()
    {
        spectrumType = RuntimeSpectrumType::RGBSpectrumType;
        nSamples = 3;
    }

    RuntimeResolutionSpectrum(Float v = 0.f)
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                runtimeSpectrum.pSS = new SampledSpectrum(v); 
                break;
            case RuntimeSpectrumType::RGBSpectrumType:
                runtimeSpectrum.pTS = new RGBSpectrum(v); 
                break;
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
        }
    }
    template <int n>
    explicit RuntimeResolutionSpectrum(const CoefficientSpectrum<n> &s)
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                runtimeSpectrum.pSS = new SampledSpectrum(s); 
                break;
            case RuntimeSpectrumType::RGBSpectrumType:
                runtimeSpectrum.pTS = new RGBSpectrum(s); 
                break;
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
        }
    }
    RuntimeResolutionSpectrum(const SampledSpectrum &s)
    {
        DCHECK(spectrumType == RuntimeSpectrumType::SampledSpectrumType);
        runtimeSpectrum.pSS = new SampledSpectrum;
        *runtimeSpectrum.pSS = s;
    }
    explicit RuntimeResolutionSpectrum(const RGBSpectrum &s)
    {
        DCHECK(spectrumType == RuntimeSpectrumType::RGBSpectrumType);
        runtimeSpectrum.pTS = new RGBSpectrum;
        *runtimeSpectrum.pTS = s;
    }

    ~RuntimeResolutionSpectrum()
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                delete runtimeSpectrum.pSS;
                break;
            case RuntimeSpectrumType::RGBSpectrumType:
                delete runtimeSpectrum.pTS;
                break;
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
        }
    }

    void Print(FILE *f) const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                runtimeSpectrum.pSS->Print(f);
                break;
            case RuntimeSpectrumType::RGBSpectrumType:
                runtimeSpectrum.pTS->Print(f);
                break;
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
        }
    }
    RuntimeResolutionSpectrum &operator+=(const RuntimeResolutionSpectrum &s2)
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                *runtimeSpectrum.pSS += *s2.runtimeSpectrum.pSS;
                return *this;
            case RuntimeSpectrumType::RGBSpectrumType:
                *runtimeSpectrum.pTS += *s2.runtimeSpectrum.pTS;
                return *this;
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return *this;
        }
    }
    RuntimeResolutionSpectrum operator+(const RuntimeResolutionSpectrum &s2) const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
            {    
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pSS = new SampledSpectrum;
                *rss.runtimeSpectrum.pSS = *runtimeSpectrum.pSS + *s2.runtimeSpectrum.pSS;
                return rss;
            }
            case RuntimeSpectrumType::RGBSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pTS = new RGBSpectrum;
                *rss.runtimeSpectrum.pTS = *runtimeSpectrum.pTS + *s2.runtimeSpectrum.pTS;
                return rss;
            }
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return *this;
        }
    }
    RuntimeResolutionSpectrum operator-(const RuntimeResolutionSpectrum &s2) const 
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pSS = new SampledSpectrum;
                *rss.runtimeSpectrum.pSS = *runtimeSpectrum.pSS - *s2.runtimeSpectrum.pSS;
                return rss;
            }
            case RuntimeSpectrumType::RGBSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pTS = new RGBSpectrum;
                *rss.runtimeSpectrum.pTS = *runtimeSpectrum.pTS - *s2.runtimeSpectrum.pTS;
                return rss;
            }       
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return *this;
        }
    }
    RuntimeResolutionSpectrum operator/(const RuntimeResolutionSpectrum &s2) const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pSS = new SampledSpectrum;
                *rss.runtimeSpectrum.pSS = *runtimeSpectrum.pSS / *s2.runtimeSpectrum.pSS;
                return rss;
            }
            case RuntimeSpectrumType::RGBSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pTS = new RGBSpectrum;
                *rss.runtimeSpectrum.pTS = *runtimeSpectrum.pTS / *s2.runtimeSpectrum.pTS;
                return rss;
            }       
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return *this;
        }
    }
    RuntimeResolutionSpectrum operator*(const RuntimeResolutionSpectrum &sp) const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pSS = new SampledSpectrum;
                *rss.runtimeSpectrum.pSS = *runtimeSpectrum.pSS * *sp.runtimeSpectrum.pSS;
                return rss;
            }
            case RuntimeSpectrumType::RGBSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pTS = new RGBSpectrum;
                *rss.runtimeSpectrum.pTS = *runtimeSpectrum.pTS * *sp.runtimeSpectrum.pTS;
                return rss;
            }  
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return *this;
        }
    }
    RuntimeResolutionSpectrum &operator*=(const RuntimeResolutionSpectrum &sp)
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                *runtimeSpectrum.pSS *= *sp.runtimeSpectrum.pSS;
                return *this;
            case RuntimeSpectrumType::RGBSpectrumType:
                *runtimeSpectrum.pTS *= *sp.runtimeSpectrum.pTS;
                return *this;
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return *this;
        }
    }
    RuntimeResolutionSpectrum operator*(Float a) const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pSS = new SampledSpectrum;
                *rss.runtimeSpectrum.pSS = *runtimeSpectrum.pSS * a;
                return rss;
            }
            case RuntimeSpectrumType::RGBSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pTS = new RGBSpectrum;
                *rss.runtimeSpectrum.pTS = *runtimeSpectrum.pTS * a;
                return rss;
            }    
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return *this;
        }
    }
    RuntimeResolutionSpectrum &operator*=(Float a)
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                *runtimeSpectrum.pSS *= a;
                return *this;
            case RuntimeSpectrumType::RGBSpectrumType:
                *runtimeSpectrum.pTS *= a;
                return *this;
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return *this;
        }
    }
    friend inline RuntimeResolutionSpectrum operator*(Float a,
                                                      const RuntimeResolutionSpectrum &s)
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pSS = new SampledSpectrum;
                *rss.runtimeSpectrum.pSS = *s.runtimeSpectrum.pSS * a;
                return rss;
            }
            case RuntimeSpectrumType::RGBSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pTS = new RGBSpectrum;
                *rss.runtimeSpectrum.pTS = *s.runtimeSpectrum.pTS * a;
                return rss;
            }
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return RuntimeResolutionSpectrum();
        }
    }
    RuntimeResolutionSpectrum operator/(Float a) const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pSS = new SampledSpectrum;
                *rss.runtimeSpectrum.pSS = *runtimeSpectrum.pSS / a;
                return rss;
            }
            case RuntimeSpectrumType::RGBSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pTS = new RGBSpectrum;
                *rss.runtimeSpectrum.pTS = *runtimeSpectrum.pTS / a;
                return rss;
            }
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return *this;
        }
    }
    RuntimeResolutionSpectrum &operator/=(Float a)
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                *runtimeSpectrum.pSS /= a;
                return *this;
            case RuntimeSpectrumType::RGBSpectrumType:
                *runtimeSpectrum.pTS /= a;
                return *this;
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return *this;
        }
    }
    bool operator==(const RuntimeResolutionSpectrum &sp) const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                return *runtimeSpectrum.pSS == *sp.runtimeSpectrum.pSS;
            case RuntimeSpectrumType::RGBSpectrumType:
                return *runtimeSpectrum.pTS == *sp.runtimeSpectrum.pTS;
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return false;
        }
    }
    bool operator!=(const RuntimeResolutionSpectrum &sp) const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                return *runtimeSpectrum.pSS != *sp.runtimeSpectrum.pSS;
            case RuntimeSpectrumType::RGBSpectrumType:
                return *runtimeSpectrum.pTS != *sp.runtimeSpectrum.pTS;
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return false;
        }
    }
    bool IsBlack() const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                return runtimeSpectrum.pSS->IsBlack();
            case RuntimeSpectrumType::RGBSpectrumType:
                return runtimeSpectrum.pTS->IsBlack();
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return false;
        }
    }
    friend RuntimeResolutionSpectrum Sqrt(const RuntimeResolutionSpectrum &s)
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pSS = new SampledSpectrum;
                *rss.runtimeSpectrum.pSS = Sqrt(*s.runtimeSpectrum.pSS);
                return rss;
            }
            case RuntimeSpectrumType::RGBSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pTS = new RGBSpectrum;
                *rss.runtimeSpectrum.pTS = Sqrt(*s.runtimeSpectrum.pTS);
                return rss;
            }
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return RuntimeResolutionSpectrum();
        }
    }
    friend inline RuntimeResolutionSpectrum Pow(const RuntimeResolutionSpectrum &s,
                                             Float e)
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pSS = new SampledSpectrum;
                *rss.runtimeSpectrum.pSS = Pow(*s.runtimeSpectrum.pSS, e);
                return rss;
            }
            case RuntimeSpectrumType::RGBSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pTS = new RGBSpectrum;
                *rss.runtimeSpectrum.pTS = Pow(*s.runtimeSpectrum.pTS, e);
                return rss;
            }            
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return RuntimeResolutionSpectrum();
        }
    }
    RuntimeResolutionSpectrum operator-() const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pSS = new SampledSpectrum;
                *rss.runtimeSpectrum.pSS = -(*runtimeSpectrum.pSS);
                return rss;
            }
            case RuntimeSpectrumType::RGBSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pTS = new RGBSpectrum;
                *rss.runtimeSpectrum.pTS = -(*runtimeSpectrum.pTS);
                return rss;
            }     
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return *this;
        }
    }
    friend RuntimeResolutionSpectrum Exp(const RuntimeResolutionSpectrum &s)
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pSS = new SampledSpectrum;
                *rss.runtimeSpectrum.pSS = Exp(*s.runtimeSpectrum.pSS);
                return rss;
            }
            case RuntimeSpectrumType::RGBSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pTS = new RGBSpectrum;
                *rss.runtimeSpectrum.pTS = Exp(*s.runtimeSpectrum.pTS);
                return rss;
            }            
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return RuntimeResolutionSpectrum();
        }
    }
    friend std::ostream &operator<<(std::ostream &os,
                                    const RuntimeResolutionSpectrum &s)
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                return os << *s.runtimeSpectrum.pSS;
            case RuntimeSpectrumType::RGBSpectrumType:
                return os << *s.runtimeSpectrum.pTS;
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return os << "ERROR";
        }
    }
    std::string ToString() const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                return runtimeSpectrum.pSS->ToString();
            case RuntimeSpectrumType::RGBSpectrumType:
                return runtimeSpectrum.pTS->ToString();
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return "ERROR";
        }
    }
    RuntimeResolutionSpectrum Clamp(Float low = 0, Float high = Infinity) const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pSS = new SampledSpectrum;
                *rss.runtimeSpectrum.pSS = runtimeSpectrum.pSS->Clamp(low, high);
                return rss;
            }
            case RuntimeSpectrumType::RGBSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pTS = new RGBSpectrum;
                *rss.runtimeSpectrum.pTS = runtimeSpectrum.pTS->Clamp(low, high);
                return rss;
            }   
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return *this;
        }
    }
    Float MaxComponentValue() const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                return runtimeSpectrum.pSS->MaxComponentValue();
            case RuntimeSpectrumType::RGBSpectrumType:
                return runtimeSpectrum.pTS->MaxComponentValue();
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return 0.0f;
        }
    }
    bool HasNaNs() const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                return runtimeSpectrum.pSS->HasNaNs();
            case RuntimeSpectrumType::RGBSpectrumType:
                return runtimeSpectrum.pTS->HasNaNs();
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return false;
        }
    }
    bool Write(FILE *f) const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                return runtimeSpectrum.pSS->Write(f);
            case RuntimeSpectrumType::RGBSpectrumType:
                return runtimeSpectrum.pTS->Write(f);
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return false;
        }
    }
    bool Read(FILE *f)
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                return runtimeSpectrum.pSS->Read(f);
            case RuntimeSpectrumType::RGBSpectrumType:
                return runtimeSpectrum.pTS->Read(f);
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return false;
        }
    }
    Float &operator[](int i)
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                return (*runtimeSpectrum.pSS)[i];
            case RuntimeSpectrumType::RGBSpectrumType:
                return (*runtimeSpectrum.pTS)[i];
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return (*runtimeSpectrum.pTS)[i];
        }
    }
    Float operator[](int i) const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                return (*runtimeSpectrum.pSS)[i];
            case RuntimeSpectrumType::RGBSpectrumType:
                return (*runtimeSpectrum.pTS)[i];
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return 0.0f;
        }
    }

    // RGBSpectrum and SampledSpectrum public functions
    RuntimeResolutionSpectrum(const RGBSpectrum &s,
                SpectrumType type = SpectrumType::Reflectance)
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                runtimeSpectrum.pSS = new SampledSpectrum(s, type); 
                break;
            case RuntimeSpectrumType::RGBSpectrumType:
                runtimeSpectrum.pTS = new RGBSpectrum(s, type); 
                break;
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
        }

    }
    static RuntimeResolutionSpectrum FromRGB(const Float rgb[3],
                          SpectrumType type = SpectrumType::Reflectance)
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
            {    
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pSS = new SampledSpectrum;
                *rss.runtimeSpectrum.pSS = SampledSpectrum::FromRGB(rgb, type);
                return rss;
            }
            case RuntimeSpectrumType::RGBSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pTS = new RGBSpectrum;
                *rss.runtimeSpectrum.pTS = RGBSpectrum::FromRGB(rgb, type);
                return rss;
            }
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return RuntimeResolutionSpectrum();
        }
    }
    void ToRGB(Float *rgb) const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                runtimeSpectrum.pSS->ToRGB(rgb);
                break;
            case RuntimeSpectrumType::RGBSpectrumType:
                runtimeSpectrum.pTS->ToRGB(rgb);
                break;
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
        }
    }
    const RGBSpectrum ToRGBSpectrum() const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                return RGBSpectrum((runtimeSpectrum.pSS)->ToRGBSpectrum());
            case RuntimeSpectrumType::RGBSpectrumType:
                return RGBSpectrum((runtimeSpectrum.pTS)->ToRGBSpectrum());
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return RGBSpectrum(0.0f);
        }
    }
    void ToXYZ(Float xyz[3]) const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                runtimeSpectrum.pSS->ToXYZ(xyz);
                break;
            case RuntimeSpectrumType::RGBSpectrumType:
                runtimeSpectrum.pTS->ToXYZ(xyz);
                break;
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
        }
    }
    static RuntimeResolutionSpectrum FromXYZ(const Float xyz[3],
                               SpectrumType type = SpectrumType::Reflectance)
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
            {    
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pSS = new SampledSpectrum;
                *rss.runtimeSpectrum.pSS = SampledSpectrum::FromXYZ(xyz, type);
                return rss;
            }
            case RuntimeSpectrumType::RGBSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pTS = new RGBSpectrum;
                *rss.runtimeSpectrum.pTS = RGBSpectrum::FromXYZ(xyz, type);
                return rss;
            }
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return RuntimeResolutionSpectrum();
        }
    }
    Float y() const
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
                return runtimeSpectrum.pSS->y();
            case RuntimeSpectrumType::RGBSpectrumType:
                return runtimeSpectrum.pTS->y();
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return 0.0f;
        }
    }
    static RuntimeResolutionSpectrum FromSampled(const Float *lambda, const Float *v, int n) 
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
            {    
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pSS = new SampledSpectrum;
                *rss.runtimeSpectrum.pSS = SampledSpectrum::FromSampled(lambda, v, n);
                return rss;
            }
            case RuntimeSpectrumType::RGBSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pTS = new RGBSpectrum;
                *rss.runtimeSpectrum.pTS = RGBSpectrum::FromSampled(lambda, v, n);
                return rss;
            }
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return RuntimeResolutionSpectrum();
        }
    }
    friend inline RuntimeResolutionSpectrum Lerp(Float t, const RuntimeResolutionSpectrum &s1, const RuntimeResolutionSpectrum &s2) 
    {
        switch(spectrumType)
        {
            case RuntimeSpectrumType::SampledSpectrumType:
            {    
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pSS = new SampledSpectrum;
                *rss.runtimeSpectrum.pSS = Lerp(t, *s1.runtimeSpectrum.pSS, *s2.runtimeSpectrum.pSS);
                return rss;
            }
            case RuntimeSpectrumType::RGBSpectrumType:
            {
                RuntimeResolutionSpectrum rss;
                rss.runtimeSpectrum.pTS = new RGBSpectrum;
                *rss.runtimeSpectrum.pTS = Lerp(t, *s1.runtimeSpectrum.pTS, *s2.runtimeSpectrum.pTS);
                return rss;
            }
            default:
                Error("Found default case in %s:%i", __FILE__, __LINE__);
                return RuntimeResolutionSpectrum();
        }
    }

    static int nSamples;
};


#endif // PBRT_RUNTIME_RESOLUTION_SPECTRUM

}  // namespace pbrt

#endif  // PBRT_CORE_SPECTRUM_H
