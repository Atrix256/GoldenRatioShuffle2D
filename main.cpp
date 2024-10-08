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

#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

typedef unsigned int uint;
typedef std::array<uint, 2> uint2;
typedef std::array<float, 2> float2;

// -1 means non deterministic
#define SEED() 435

#define NUM_FRAMES() 10  // 0 to save all frames
static const uint2 c_imageSizes2D[] = { { 64, 64 }, { 128, 128 }, { 256, 256 }, { 512, 512 }, { 1024, 1024 }, { 2048, 2048 } };

static const uint c_numIntegrationtests = 1000;

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

float Lerp(float A, float B, float t)
{
	return A * (1.0f - t) + B * t;
}

float Frac(float v)
{
	return v - std::floor(v);
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
			// Save the last frame
			saveFrame |= (index == dims[0] * dims[1] - 1); 

			// save frames periodically
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

	// 1D Shuffle
	// Use a standard 1D LDS Shuffle, and just convert that to 2d (index % width, index / width)
	{
		LDShuffle shuffle(dims[0] * dims[1], seed);
		DoTest2D_SingleTest("1DShuffle", dims,
			[&shuffle, &dims](uint index)
			{
				// Get the shuffled index
				uint shuffledIndex = shuffle.GetValueAtIndex(index);

				// Convert that to a 2d coordinate
				uint x = shuffledIndex % dims[0];
				uint y = shuffledIndex / dims[0];

				return uint2{ x, y };
			}
		);
	}

	// White noise
	{
		std::mt19937 rng(seed);

		std::vector<uint2> points(dims[0] * dims[1]);
		for (uint i = 0; i < (uint)points.size(); ++i)
			points[i] = uint2{ i % dims[0], i / dims[0] };
		std::shuffle(points.begin(), points.end(), rng);

		DoTest2D_SingleTest("White", dims,
			[&points](uint index)
			{
				return points[index];
			}
		);
	}

	/*
	// Random (white noise, but not a shuffle!)
	{
		std::mt19937 rng(seed);

		std::uniform_int_distribution<int> dist(0, dims[0] - 1);

		DoTest2D_SingleTest("Rand", dims,
			[&dist, &rng](uint index)
			{
				return uint2{ (uint)dist(rng),(uint)dist(rng) };
			}
		);
	}
	*/

	/*
	// R2 sequence - V1
	// Do the normal R2 sequence, but convert it to integers
	{
		std::mt19937 rng(seed);
		std::uniform_real_distribution<float> dist(0.0f, 1.0f);

		float2 point = { dist(rng), dist(rng) };

		DoTest2D_SingleTest("R2v1", dims,
			[&point, &dims](uint index)
			{
				uint2 pointi;
				pointi[0] = Clamp<uint>(uint(point[0] * float(dims[0]) + 0.5f), 0, dims[0] - 1);
				pointi[1] = Clamp<uint>(uint(point[1] * float(dims[1]) + 0.5f), 0, dims[1] - 1);

				point[0] = Frac(point[0] + c_R2_a1);
				point[1] = Frac(point[1] + c_R2_a2);

				return pointi;
			}
		);
	}
	*/

	/*
	// R2 sequence - V2
	// Run a shuffle for each axis, using the R2 sequence
	// Total failure - the shuffles repeat every dims[].
	{
		LDShuffle shuffle1(dims[0], seed, c_R2_a1);
		LDShuffle shuffle2(dims[1], seed, c_R2_a2);

		DoTest2D_SingleTest("R2v2", dims,
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
	// R2 sequence - V3
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

		DoTest2D_SingleTest("R2v3", dims,
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
	// Now, repeat that, but with an offset to X, using another 1D LDS shuffle over width.
	// It gives the same point set, but being offset it was hard to tell and was decently spaced.
	// I also had to make sure the various numbers were coprime to each other, and the dimensions of the texture.
	// I also had to use 3 irrationals: golden ratio, sqrt(2), sqrt(3).
	// This worked pretty decently at some screen resolutions but not others.
}

void DoIntegrationTest(uint seed)
{
	int w, h, c;
	unsigned char* pixels = stbi_load("cabin.png", &w, &h, &c, 1);

	// calculate the actual average value of the image
	float average = 0.0f;
	for (size_t i = 0; i < 512 * 512; ++i)
		average = Lerp(average, float(pixels[i]) / 255.0f, 1.0f / float(i + 1));
	//printf("Actual average: %f\n", average);

	printf("Integrating....\n");

	std::mt19937 rng(seed);

	std::vector<float> mse[4];
	mse[0].resize(512 * 512, 0.0f);
	mse[1].resize(512 * 512, 0.0f);
	mse[2].resize(512 * 512, 0.0f);
	mse[3].resize(512 * 512, 0.0f);

	int lastPercent = -1;
	for (uint testIndex = 0; testIndex < c_numIntegrationtests; ++testIndex)
	{
		int percent = int(100.0f * float(testIndex) / float(c_numIntegrationtests - 1));
		if (percent != lastPercent)
		{
			lastPercent = percent;
			printf("\r%i%%", percent);
		}

		// Make a shuffled list of pixels to get a white noise ordering
		std::vector<uint> shuffleOrder(512 * 512);
		{
			for (uint i = 0; i < shuffleOrder.size(); ++i)
				shuffleOrder[i] = i;
			std::shuffle(shuffleOrder.begin(), shuffleOrder.end(), rng);
		}

		// One low discrepancy shuffler will be used by both the Hilbert and ZOrder versions
		LDShuffle shuffle(512 * 512, rng());

		float avgs[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
		for (size_t sampleIndex = 0; sampleIndex < 512 * 512; ++sampleIndex)
		{
			// white noise
			avgs[0] = Lerp(avgs[0], float(pixels[shuffleOrder[sampleIndex]]) / 255.0f, 1.0f / float(sampleIndex + 1));

			// Get the shuffled index
			uint shuffledIndex = shuffle.GetValueAtIndex((uint)sampleIndex);

			// Hilbert
			{
				int x, y;
				Hilbert::d2xy(512 * 512, shuffledIndex, &x, &y);

				uint readIndex = y * 512 + x;
				avgs[1] = Lerp(avgs[1], float(pixels[readIndex]) / 255.0f, 1.0f / float(sampleIndex + 1));
			}

			// ZOrder
			{
				uint x, y;
				ZOrder::OneDToTwoD(shuffledIndex, x, y);

				uint readIndex = y * 512 + x;
				avgs[2] = Lerp(avgs[2], float(pixels[readIndex]) / 255.0f, 1.0f / float(sampleIndex + 1));
			}

			// 1D Shuffler without a curve
			{
				uint readIndex = shuffledIndex;
				avgs[3] = Lerp(avgs[3], float(pixels[readIndex]) / 255.0f, 1.0f / float(sampleIndex + 1));
			}

			// Keep track of MSE across tests
			for (int noiseIndex = 0; noiseIndex < 4; ++noiseIndex)
			{
				float error = std::abs(avgs[noiseIndex] - average);
				mse[noiseIndex][sampleIndex] = Lerp(mse[noiseIndex][sampleIndex], error * error, 1.0f / float(testIndex + 1));
			}
		}
	}
	printf("\r100%%");

	// Write the data out
	FILE* file = nullptr;
	fopen_s(&file, "out/_integration.csv", "wb");
	fprintf(file, "\"Index\",\"White\",\"Hilbert\",\"ZOrder\",\"1DShuffler\"\n");
	for (uint sampleIndex = 0; sampleIndex < 512 * 512; sampleIndex += 512)
	{
		fprintf(file, "\"%u\",\"%f\",\"%f\",\"%f\",\"%f\"\n", sampleIndex+1, std::sqrt(mse[0][sampleIndex]), std::sqrt(mse[1][sampleIndex]), std::sqrt(mse[2][sampleIndex]), std::sqrt(mse[3][sampleIndex]));
	}
	fclose(file);


	stbi_image_free(pixels);
}

int main(int argc, char** argv)
{
	_mkdir("out");

	// initialize RNG
	std::random_device rd;
	uint seed = (SEED() == -1) ? rd() : SEED();
	printf("Seed = %u\n\n", seed);

	for (const uint2& size : c_imageSizes2D)
	{
		printf("============= 2D Tests %u x %u =============\n\n", size[0], size[1]);
		DoTest2D(size, seed);
	}

	DoIntegrationTest(seed);

	return 0;
}
