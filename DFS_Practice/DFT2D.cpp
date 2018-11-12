// DFT_2D (Image Processing)

#define CRTDBG_MAP_ALLOC

#include "DFT2D.h"

#define HEIGHT 256
#define WIDTH 256

#define BUTTERWORTH_ORDER 2

typedef unsigned char BYTE;

template<typename T> T** MemAlloc2D(int nHeight, int nWidth, int nInitVal)
{
	T** rtn = new T*[nHeight];
	for (int h = 0; h < nHeight; h++)
	{
		rtn[h] = new T[nWidth];
		memset(rtn[h], nInitVal, sizeof(T) * nWidth);
	}
	return rtn;
}
template<typename T> void MemFree2D(T** arr2D, int nHeight)
{
	for (int h = 0; h < nHeight; h++)
	{
		delete[] arr2D[h];
	}
	delete[] arr2D;
}
void FileRead(const char* strFilename, BYTE** arr2D, int nHeight, int nWidth)
{
	FILE* fp_in = fopen(strFilename, "rb");
	for (int h = 0; h < nHeight; h++)
	{
		fread(arr2D[h], sizeof(BYTE), nWidth, fp_in);
	}

	fclose(fp_in);
}
void FileWrite(const char* strFilename, BYTE** arr2D, int nHeight, int nWidth)
{
	FILE* fp_out = fopen(strFilename, "wb");
	for (int h = 0; h < nHeight; h++)
	{
		fwrite(arr2D[h], sizeof(BYTE), nWidth, fp_out);
	}

	fclose(fp_out);
}
BYTE clip_d(double a)
{
	if (a > 255)
	{
		return (BYTE)255;
	}
	if (a < 0)
	{
		return 0;
	}
	else return (BYTE)floor(a + 0.5);
}

void DFT_2D(double **dReal_2D, double **dImag_2D, int nHeight, int nWidth, int DFT_ID)
{
	double** tempReal_2D = MemAlloc2D<double>(HEIGHT, WIDTH, 0);
	double** tempImag_2D = MemAlloc2D<double>(HEIGHT, WIDTH, 0);
	for (int h = 0; h < nHeight; h++)
	{
		DFT_1D(dReal_2D[h], dImag_2D[h], nWidth, DFT_ID);
	}
	for (int h = 0; h < nHeight; h++)
	{
		for (int w = 0; w < nWidth; w++)
		{
			tempReal_2D[h][w] = dReal_2D[w][h];
			tempImag_2D[h][w] = dImag_2D[w][h];
		}
	}
	for (int h = 0; h < nHeight; h++)
	{
		memcpy(dReal_2D[h], tempReal_2D[h], sizeof(double) * nWidth);
		memcpy(dImag_2D[h], tempImag_2D[h], sizeof(double) * nWidth);
	}
	for (int h = 0; h < nHeight; h++)
	{
		DFT_1D(dReal_2D[h], dImag_2D[h], nWidth, DFT_ID);
	}
}

void DFT_1D(double* dReal_1D, double* dImag_1D, int nLength, int DFT_ID)
{
	double dArg, dCosArg, dSinArg;

	double *dTemp_R = new double[nLength];
	double *dTemp_I = new double[nLength];

	memset(dTemp_R, 0, sizeof(double) * nLength);
	memset(dTemp_I, 0, sizeof(double) * nLength);

	for (int m = 0; m < nLength; m++)
	{
		dTemp_R[m] = 0;
		dTemp_I[m] = 0;
		dArg = DFT_ID * 2.0 * PI * (double)m / (double)nLength;
		for (int n = 0; n < nLength; n++)
		{
			dCosArg = cos(n * dArg);
			dSinArg = sin(n * dArg);

			dTemp_R[m] += (dReal_1D[n] * dCosArg - dImag_1D[n] * dSinArg);
			dTemp_I[m] += (dReal_1D[n] * dSinArg + dImag_1D[n] * dCosArg);
		}
	}

	if (DFT_ID == IDFT)
	{
		for (int m = 0; m < nLength; m++)
		{
			dReal_1D[m] = dTemp_R[m] / (double)nLength;
			dImag_1D[m] = dTemp_I[m] / (double)nLength;
		}
	}
	else
	{
		for (int m = 0; m < nLength; m++)
		{
			dReal_1D[m] = dTemp_R[m];
			dImag_1D[m] = dTemp_I[m];
		}
	}

	delete[] dTemp_R;
	delete[] dTemp_I;

}
void ConvImage(BYTE** arr2D, double** dReal_2D, double** dImag_2D, int nHeight, int nWidth, int DFT_ID)
{
	double dMax = 0, dNor, dMean = 0;
	double** dTemp_2D = MemAlloc2D<double>(nHeight, nWidth, 0);

	for (int h = 0; h < nHeight; h++)
	{
		memcpy(dTemp_2D[h], dReal_2D[h], sizeof(double) * nWidth);
	}

	if (DFT_ID == DFT)
	{
		for (int h = 0; h < nHeight; h++)
		{
			for (int w = 0; w < nWidth; w++)
			{
				double dTemp = sqrt(dReal_2D[h][w] * dReal_2D[h][w] + dImag_2D[h][w] * dImag_2D[h][w]);
				dTemp_2D[h][w] = log10(dTemp + 1);

				if (dMax < dTemp_2D[h][w])
				{
					dMax = dTemp_2D[h][w];
				}
			}
		}

		dNor = 255 / dMax;
	}
	else
	{
		dNor = 1;
	}

	for (int h = 0; h < nHeight; h++)
	{
		for (int w = 0; w < nWidth; w++)
		{
			double dTemp = (dTemp_2D[h][w] * dNor);
			arr2D[h][w] = clip_d(dTemp);
		}
	}
	if (DFT_ID == DFT)
	{
		SetCenterDC<BYTE>(arr2D, nHeight, nWidth);
	}

	MemFree2D<double>(dTemp_2D, nHeight);
}
template <typename T> void SetCenterDC(T** arr2D, int nHeight, int nWidth)
{
	int nHalf_H = nHeight / 2, nHalf_W = nWidth / 2;

	T** arr2DTemp = MemAlloc2D<T>(nHeight, nWidth, 0);

	for (int h = 0; h < nHeight; h++)
	{
		for (int w = 0; w < nWidth; w++)
		{
			arr2DTemp[h][w] = arr2D[h][w];
		}
	}

	for (int h = 0; h < nHalf_H; h++)
	{
		for (int w = 0; w < nHalf_W; w++)
		{
			arr2D[h][w] = arr2DTemp[h + nHalf_H][w + nHalf_W];
			arr2D[h + nHalf_H][w] = arr2DTemp[h][w + nHalf_W];
			arr2D[h][w + nHalf_W] = arr2DTemp[h + nHalf_H][w];
			arr2D[h + nHalf_H][w + nHalf_W] = arr2DTemp[h][w];
		}
	}

	MemFree2D<T>(arr2DTemp, nHeight);
}
void LowPassIdeal(int nHeight, int nWidth, int nThres, double** dFilter)
{
	int w2 = nWidth / 2;
	int h2 = nHeight / 2;
	int temp_h, temp_w;

	for (int j = 0; j < nHeight; j++)
	{
		for (int i = 0; i < nWidth; i++)
		{
			temp_w = i + w2;
			temp_h = j + h2;

			if (temp_w >= w2) temp_w -= w2;
			if (temp_h >= h2) temp_h -= h2;


			double check = sqrt((double)((temp_h - h2)*(temp_h - h2)) + ((temp_w - w2)*(temp_w - w2)));

			if (check > nThres)
			{
				dFilter[j][i] = 0;
			}
			else
				dFilter[j][i] = 1;

		}
	}
}

void LowPassGaussian(int nHeight, int nWidth, int nThres, double** dFilter)
{
	int w2 = nWidth / 2;
	int h2 = nHeight / 2;
	int temp_h, temp_w;

	double dist2, hval;

	for (int j = 0; j < nHeight; j++)
	{
		for (int i = 0; i < nWidth; i++)
		{
			temp_w = i + w2;
			temp_h = j + h2;

			if (temp_w >= w2) temp_w -= w2;
			if (temp_h >= h2) temp_h -= h2;


			dist2 = (double)((temp_h - h2)*(temp_h - h2)) + ((temp_w - w2)*(temp_w - w2));

			hval = exp(-dist2 / (2 * nThres*nThres));

			dFilter[j][i] = hval;

		}
	}
}

void LowPassButterworth(int nHeight, int nWidth, int nThres, double** dFilter)
{
	int w2 = nWidth / 2;
	int h2 = nHeight / 2;
	int temp_h, temp_w;

	double hval;
	int n = 1; //order

	for (int j = 0; j < nHeight; j++)
	{
		for (int i = 0; i < nWidth; i++)
		{
			temp_w = i + w2;
			temp_h = j + h2;

			if (temp_w >= w2) temp_w -= w2;
			if (temp_h >= h2) temp_h -= h2;


			double check = sqrt((double)((temp_h - h2)*(temp_h - h2)) + ((temp_w - w2)*(temp_w - w2)));

			hval = 1.0 + pow((check / nThres), 2 * n);

			hval = 1 / hval;

			dFilter[j][i] = hval;

		}
	}
}
void LowPassFilter(double** dReal2D, double** dImag2D, int nHeight, int nWidth, int nThres, int nFilterType)
{
	int nHalf_H = nHeight / 2, nHalf_W = nWidth / 2;

	double** dFilter = MemAlloc2D<double>(nHeight, nWidth, 0);

	SetCenterDC<double>(dReal2D, nHeight, nWidth);
	SetCenterDC<double>(dImag2D, nHeight, nWidth);

	// Make Filter
	switch (nFilterType)
	{
	case LPF_IDEAL:
		LowPassIdeal(nHeight, nWidth, nThres, dFilter);

		break;
	case LPF_BUTTERWORTH:
		LowPassButterworth(nHeight, nWidth, nThres, dFilter);

		break;
	case LPF_GAUSSIAN:
		LowPassGaussian(nHeight, nWidth, nThres, dFilter);

		break;
	}

	for (int h = 0; h < nHeight; h++)
	{
		for (int w = 0; w < nWidth; w++)
		{
			dReal2D[h][w] = dReal2D[h][w] * dFilter[h][w];
			dImag2D[h][w] = dImag2D[h][w] * dFilter[h][w];
		}
	}
	SetCenterDC<double>(dReal2D, nHeight, nWidth);
	SetCenterDC<double>(dImag2D, nHeight, nWidth);

	MemFree2D<double>(dFilter, nHeight);
}
void HighPassFilter(double** dReal2D, double** dImag2D, int nHeight, int nWidth, int nThres, int nFilterType)
{
	int nHalf_H = nHeight / 2, nHalf_W = nWidth / 2;

	double** dFilter = MemAlloc2D<double>(nHeight, nWidth, 0);

	SetCenterDC<double>(dReal2D, nHeight, nWidth);
	SetCenterDC<double>(dImag2D, nHeight, nWidth);

	// Make Filter
	switch (nFilterType)
	{
	case HPF_IDEAL:

		break;
	case HPF_BUTTERWORTH:

		break;
	case HPF_GAUSSIAN:

		break;
	}

	for (int h = 0; h < nHeight; h++)
	{
		for (int w = 0; w < nWidth; w++)
		{
			dReal2D[h][w] = dReal2D[h][w] * dFilter[h][w];
			dImag2D[h][w] = dImag2D[h][w] * dFilter[h][w];
		}
	}
	SetCenterDC<double>(dReal2D, nHeight, nWidth);
	SetCenterDC<double>(dImag2D, nHeight, nWidth);
}