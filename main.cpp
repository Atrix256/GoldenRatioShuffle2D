#define _CRT_SECURE_NO_WARNINGS // for stb

#include <stdio.h>
#include <cmath>
#include <random>
#include <vector>
#include <array>
#include <direct.h>
#include "Hilbert.h"
#include "LDShuffle.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

typedef unsigned int uint;
typedef std::array<uint, 2> uint2;

// -1 means non deterministic
#define SEED() 1

#define PRINT_FRAMES_2D() true
#define NUM_FRAMES_2D() 30  // 0 to save all frames
static const uint2 c_imageSizes2D[] = { { 64, 64 } };//, { 128, 128 }, { 256, 256 }, { 512, 512 }, { 640, 480 } };

template <typename T>
T Frac(T f)
{
	return std::fmod(f, T(1.0));
}

template <typename T>
T Clamp(T value, T themin, T themax)
{
	if (value <= themin)
		return themin;
	else if (value >= themax)
		return themax;
	else
		return value;
}

template <typename LAMBDA>
void DoTest2D_Single(const uint2& dims, const LAMBDA& lambda)
{
	std::vector<int> visitCount(dims[0]*dims[1], 0);

	int errorCount = 0;
	for (uint index = 0; index < dims[0] * dims[1]; ++index)
	{
		float percent = float(index) / float(dims[0] * dims[1] - 1);
		uint2 shuffleItem = lambda(index, percent);

		uint shuffleItemFlat = shuffleItem[1] * dims[0] + shuffleItem[0];

		visitCount[shuffleItemFlat]++;
		if (visitCount[shuffleItemFlat] > 1)
			errorCount++;
	}

	printf("%i Duplicates (%0.2f%%)\n\n", errorCount, 100.0f * float(errorCount) / float(dims[0] * dims[1]));
}

void DoTest2D(const uint2& dims, uint seed)
{
	printf("========== 2D: %u x %u ==========\n\n", dims[0], dims[1]);

	// TODO: make this work. only save 30 frames, for example. it'll make a nice gif or something.
	// TODO: maybe get rid of video stuff, and just make gifs?
	uint numVideoFrames = ((NUM_FRAMES_2D() > 0) ? NUM_FRAMES_2D() : dims[0] * dims[1]);

	// Hilbert
	std::vector<unsigned char> image(dims[0] * dims[1], 0);
	std::vector<unsigned char> image2(dims[0] * dims[1], 0);

	LDShuffle shuffle(dims[0] * dims[1], seed);

	printf("Hilbert (%u): ", shuffle.GetCoprime());

	DoTest2D_Single(dims,
		[&dims, &image, &image2, &shuffle](uint index, float percent)
		{
			// Get the shuffled index
			uint shuffledIndex = shuffle.GetValueAtIndex(index);

			// Convert that to a 2d coordinate using the hilbert curve
			int x = 0;
			int y = 0;
			Hilbert::d2xy(dims[0] * dims[1], shuffledIndex, &x, &y);

			// make the return value
			uint2 ret = { uint(x), uint(y) };

			// draw the point
			image[ret[1] * dims[0] + ret[0]] = (unsigned char)Clamp(percent * 255.0f, 0.0f, 255.0f);
			image2[ret[1] * dims[0] + ret[0]] = 255;

			if (PRINT_FRAMES_2D())
			{
				char fileName[256];
				sprintf_s(fileName, "out/Hilbert_%u_%u_%u.png", dims[0], dims[1], index + 1);
				stbi_write_png(fileName, (int)dims[0], (int)dims[1], 1, image2.data(), 0);
			}

			return ret;
		}
	);

	char fileName[256];
	sprintf_s(fileName, "out/Hilbert_%u_%u.png", dims[0], dims[1]);
	stbi_write_png(fileName, (int)dims[0], (int)dims[1], 1, image.data(), 0);
}

int main(int argc, char** argv)
{
	_mkdir("out");

	// initialize RNG
	std::random_device rd;
	uint seed = (SEED() == -1) ? rd() : SEED();
	printf("Seed = %u\n\n", seed);

	for (const uint2& size : c_imageSizes2D)
		DoTest2D(size, seed);

	return 0;
}


/*
TODO:
- need to clean all this stuff up
- should we integrate an image and graph white noise (shuffled) vs this? 
- show how to invert it
- write blog post
 - include DFT (from gigi?)

*/
