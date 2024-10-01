#define _CRT_SECURE_NO_WARNINGS // for stb

#include <stdio.h>
#include <random>
#include <vector>
#include <array>
#include <direct.h>
#include "Hilbert.h"
#include "ZOrder.h"
#include "LDShuffle.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

typedef unsigned int uint;
typedef std::array<uint, 2> uint2;

// -1 means non deterministic
#define SEED() 435

#define NUM_FRAMES() 10  // 0 to save all frames
static const uint2 c_imageSizes2D[] = { { 64, 64 }, { 128, 128 }, { 256, 256 }, { 512, 512 }, { 1024, 1024 }, { 2048, 2048 } };

// R2 is from http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
static const float c_R2_g = 1.32471795724474602596f;
static const float c_R2_a1 = 1 / c_R2_g;
static const float c_R2_a2 = 1 / (c_R2_g * c_R2_g);

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

// From https://blog.demofox.org/2015/01/24/programmatically-calculating-gcd-and-lcm/
static inline uint CalculateGCD(uint smaller, uint larger)
{
	// make sure A <= B before starting
	if (larger < smaller)
		std::swap(smaller, larger);

	// loop
	while (1)
	{
		// if the remainder of larger / smaller is 0, they are the same
		// so return smaller as the GCD
		uint remainder = larger % smaller;
		if (remainder == 0)
			return smaller;

		// otherwise, the new larger number is the old smaller number, and
		// the new smaller number is the remainder
		larger = smaller;
		smaller = remainder;
	}
}

// Returns true, if c is coprime to all of the values
bool IsCoprime(uint c, const std::vector<uint>& values)
{
	for (uint value : values)
	{
		if (CalculateGCD(c, value) != 1)
			return false;
	}
	return true;
}

// Find the coprime of <values> that is nearest to target
uint GetNearestCoprime(uint target, const std::vector<uint>& values)
{
	uint coprime = 0;

	uint offset = 0;
	while (1)
	{
		if (offset < target)
		{
			coprime = target - offset;
			if (IsCoprime(coprime, values))
				break;
		}

		coprime = target + offset + 1;
		if (IsCoprime(coprime, values))
			break;

		offset++;
	}

	return coprime;
}

template <typename LAMBDA>
void DoTest2D_SingleTest(const char* testName, const uint2& dims, const LAMBDA& lambda)
{
	std::vector<unsigned char> image(dims[0] * dims[1], 0);
	std::vector<unsigned char> image2(dims[0] * dims[1], 0);

	std::vector<int> visitCount(dims[0] * dims[1], 0);

	int errorCount = 0;
	int frameIndex = 0;

	printf("%s %u x %u: ", testName, dims[0], dims[1]);

	for (uint index = 0; index < dims[0] * dims[1]; ++index)
	{
		// Get the point
		uint2 p = lambda(index);

		// draw the point
		float percent = float(index) / float(dims[0] * dims[1] - 1);
		image[p[1] * dims[0] + p[0]] = (unsigned char)Clamp(percent * 255.0f, 0.0f, 255.0f);
		image2[p[1] * dims[0] + p[0]] = 255;

		// determine whether we should save this frame out
		bool saveFrame = false;
		if (NUM_FRAMES() == 0)
		{
			saveFrame = true;
		}
		else
		{
			saveFrame |= (index == 0); // Save the first frame
			saveFrame |= (index == dims[0] * dims[1] - 1); // Save the last frame

			if (index > 0)
			{
				int lastSection = (index - 1) * NUM_FRAMES() / (dims[0] * dims[1]);
				int thisSection = index * NUM_FRAMES() / (dims[0] * dims[1]);
				saveFrame |= (lastSection != thisSection);
			}
		}

		// save it if so
		if (saveFrame)
		{
			char fileName[256];
			sprintf_s(fileName, "out/%s_%u_%u_%u.png", testName, dims[0], dims[1], frameIndex);
			stbi_write_png(fileName, (int)dims[0], (int)dims[1], 1, image2.data(), 0);
			frameIndex++;
		}

		// Check for hitting any values more than once
		uint shuffleItemFlat = p[1] * dims[0] + p[0];
		visitCount[shuffleItemFlat]++;
		if (visitCount[shuffleItemFlat] > 1)
			errorCount++;
	}

	// Write out the ordering of the texture
	char fileName[256];
	sprintf_s(fileName, "out/%s_%u_%u.png", testName, dims[0], dims[1]);
	stbi_write_png(fileName, (int)dims[0], (int)dims[1], 1, image.data(), 0);

	// report how many duplicates were encountered
	printf("%i Duplicates (%0.2f%%)\n\n", errorCount, 100.0f * float(errorCount) / float(dims[0] * dims[1]));
}

void DoTest2D(const uint2& dims, uint seed)
{
	// Hilbert
	// Use a standard 1D LDS Shuffle, but convert the 1d index to 2d using a hilbert curve.
	{
		LDShuffle shuffle(dims[0] * dims[1], seed);
		DoTest2D_SingleTest("Hilbert", dims,
			[&shuffle, &dims] (uint index)
			{
				// Get the shuffled index
				uint shuffledIndex = shuffle.GetValueAtIndex(index);

				// Convert that to a 2d coordinate using the hilbert curve
				int x = 0;
				int y = 0;
				Hilbert::d2xy(dims[0] * dims[1], shuffledIndex, &x, &y);

				// Verify that we can convert it back to a 1d value and that it's the same value
				uint shuffledIndexRoundTrip = Hilbert::xy2d(dims[0] * dims[1], x, y);
				if (shuffledIndexRoundTrip != shuffledIndex)
					printf("Hilbert inversion failed!\n");

				return uint2{ (uint)x, (uint)y };
			}
		);
	}

	// Z order
	// Use a standard 1D LDS Shuffle, but convert the 1d index to 2d using a Z order curve.
	{
		LDShuffle shuffle(dims[0] * dims[1], seed);
		DoTest2D_SingleTest("ZOrder", dims,
			[&shuffle](uint index)
			{
				// Get the shuffled index
				uint shuffledIndex = shuffle.GetValueAtIndex(index);

				// Convert that to a 2d coordinate using the Z order curve
				uint x = 0;
				uint y = 0;
				ZOrder::OneDToTwoD(shuffledIndex, x, y);

				// Verify that we can convert it back to a 1d value and that it's the same value
				uint shuffledIndexRoundTrip = ZOrder::TwoDToOneD(x, y);
				if (shuffledIndexRoundTrip != shuffledIndex)
					printf("ZOrder inversion failed!\n");

				return uint2{ x, y };
			}
		);
	}

	/*
	// R2 sequence - Naive
	// Run a shuffle for each axis, using the R2 sequence
	// Total failure - the shuffles repeat every dims[].
	{
		LDShuffle shuffle1(dims[0], seed, c_R2_a1);
		LDShuffle shuffle2(dims[1], seed, c_R2_a2);

		DoTest2D_SingleTest("R2Naive", dims,
			[&shuffle1, &shuffle2](uint index)
			{
				uint x = shuffle1.GetValueAtIndex(index);
				uint y = shuffle2.GetValueAtIndex(index);
				return uint2{ x, y };
			}
		);
	}
	*/

	/*
	// R2 sequence - Better
	// Convert the 2d R2 constants into a 1D offset index.
	// Run a 1D shuffle based on that, and use that 1D index as the 2d point on the image without transformation.
	// Very streaky unfortunately.
	{
		uint2 R2Offset = {
			uint(c_R2_a1 * float(dims[0]) + 0.5f),
			uint(c_R2_a2 * float(dims[1]) + 0.5f)
		};

		uint R2OffsetFlat = R2Offset[1] * dims[0] + R2Offset[0];

		uint coprime = GetNearestCoprime(R2OffsetFlat, { dims[0] * dims[1] });

		LDShuffle shuffle(dims[0] * dims[1], seed, coprime);

		DoTest2D_SingleTest("R2Better", dims,
			[&shuffle, &dims](uint index)
			{
				uint shuffledIndex = shuffle.GetValueAtIndex(index);
				uint x = shuffledIndex % dims[0];
				uint y = shuffledIndex / dims[0];
				return uint2{ x, y };
			}
		);
	}
	*/

	// I tried some other things, but nothing really worked out.
	// One thing I tried was to do a 1D shuffle for width to get all X values exactly once.
	// While doing that, I also did a 1D shuffle for height to get all the Y values exactly once.
	// At the end of this, there was one pixel per row.
	// Now, repeat that, but with an offset to X, using another 1D LDS shuffle.
	// I also had to make sure the various numbers were coprime to each other, and the dimensions of the texture.
	// I also had to use 3 irrationals: golden ratio, sqrt(2), sqrt(3).
	// This worked pretty decently at some screen resolutions but not others.
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
- could also add r2 back in.
- could also try z order. You just interleave bits: https://en.wikipedia.org/wiki/Z-order_curve
 - can re-centralize tests. make it easier for people to test their own?
- why does hilbert 1024 look so strange, with clumps? 4096 looks ok. 2048 looks strange too. maybe just note it for now and move on.
 - hilbert is such that nearby 1d points are nearby in 2d. The converse may not be true, that distance 1d points are always distant in 2D?
- should we integrate an image and graph white noise (shuffled) vs this?  maybe with multiple seeds
- show how to invert it.
- write blog post
 - include DFT (from gigi? or python snippets. whatever)
 - needs to be a power of 2, how you wrote it. would need to adapt the hilbert code for non power of 2
 - or break the image into powers of 2 and shuffle each individually, but use a weighted round robin to visit each? (could link to low discrepancy weighted round robin)

*/
