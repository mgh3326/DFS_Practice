#pragma once

#include <crtdbg.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.1415926535
#define MAX(a,b)		(((a)>(b)) ? (a) : (b))
#define DFT -1
#define IDFT 1

typedef unsigned char BYTE;

enum FilterType
{
	LPF_IDEAL,
	LPF_BUTTERWORTH,
	LPF_GAUSSIAN,
	HPF_IDEAL,
	HPF_BUTTERWORTH,
	HPF_GAUSSIAN
};


// Memory Management
template<typename T> T** MemAlloc2D(int nHeight, int nWidth, int nInitVal);
template<typename T> void MemFree2D(T** arr2D, int nHeight);

// File Read & Write
void FileRead(const char* strFilename, BYTE** arr2D, int nHeight, int nWidth);
void FileWrite(const char* strFilename, BYTE** arr2D, int nHeight, int nWidth);

// Clipping (double)
BYTE clip_d(double a);

// 2-D DFT & Inverse 2-D DFT
void DFT_2D(double **dReal_2D, double **dImag_2D, int nHeight, int nWidth, int DFT_ID);
void DFT_1D(double* dReal_1D, double* dImag_1D, int nLength, int DFT_ID);

// Visualization of 2-D DFT Magnitude
void ConvImage(BYTE** arr2D, double** dReal_2D, double** dImag_2D, int nHeight, int nWidth, int DFT_ID);
template <typename T> void SetCenterDC(T** arr2D, int nHeight, int nWidth);

// Filtering In Freq Domain
void LowPassFilter(double** dReal2D, double** dImag2D, int nHeight, int nWidth, int nRadius, int nFilterType);
void HighPassFilter(double** dReal2D, double** dImag2D, int nHeight, int nWidth, int nThres, int nFilterType);
