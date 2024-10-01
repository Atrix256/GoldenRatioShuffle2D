#pragma once

// Z order curve, or Morton order curve
// A 1D integer becomes 2D by using the even bits for X and odd bits for Y
// a 2D integer vector becomes 1D by interleaving the bits

namespace ZOrder
{
	inline void OneDToTwoD(unsigned int index, unsigned int& x, unsigned int& y)
	{
		unsigned int shift = 0;
		x = y = 0;
		while (index != 0)
		{
			x |= (index & 1) << shift;
			index = index >> 1;

			y |= (index & 1) << shift;
			index = index >> 1;

			shift++;
		}
	}

	inline unsigned int TwoDToOneD(unsigned int x, unsigned int y)
	{
		unsigned int ret = 0;
		unsigned int shift = 0;
		while (x != 0 || y != 0)
		{
			ret |= ((x & 1) << shift);
			x = x >> 1;
			shift++;

			ret |= ((y & 1) << shift);
			y = y >> 1;
			shift++;
		}

		return ret;
	}
};
