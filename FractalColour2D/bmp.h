#pragma once
#include "stdafx.h"

BYTE* ConvertRGBToBMPBuffer(BYTE* Buffer, int width, int height, long* newsize);
bool SaveBMP(BYTE* Buffer, int width, int height, long paddedsize, LPCTSTR bmpfile);
BYTE* LoadBMP(int* width, int* height, long* size, LPCTSTR bmpfile);

