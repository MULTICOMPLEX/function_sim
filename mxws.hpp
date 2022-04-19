
#include <random>
#include <concepts>

class mxws
{

private:

	std::random_device r;

public:

	uint64_t x, w;

	typedef uint32_t result_type;

	void seed()
	{
		init();
	}

	mxws(const std::seed_seq& seq)
	{
		if (seq.size() == 2)
		{
			std::vector<uint32_t> seeds(seq.size());
			seq.param(seeds.rbegin());
			w = (uint64_t(seeds[1]) << 32) | seeds[0];
			x = 1;
		}

		else init();
	}

	mxws()
	{
		w = (uint64_t(r()) << 32) | r();
		x = 1;
	}

	mxws(const uint64_t& seed)
	{
		init();
	}

	void init()
	{
		w = (uint64_t(r()) << 32) | r();
		x = 1;
	}

	void init(const uint64_t& seed)
	{
		w = seed;
		x = 1;
	}

	virtual ~mxws() = default;

	static uint32_t min() { return std::numeric_limits<uint32_t>::min(); }
	static uint32_t max() { return std::numeric_limits<uint32_t>::max(); }

	inline uint32_t operator()()
	{
		x *= w;
		x = (x >> 32) | (x << 32);
		w += x;
		return uint32_t(x);
	}

	template <typename T>
		requires std::floating_point<T>
	inline T operator()(const T& f)
	{
		return (*this)() / 4294967296. * f;
	}

	template <typename T>
		requires std::floating_point<T>
	inline T operator()(const T& min, const T& max)
	{
		return (*this)(1.0) * (max - min) + min;
	}

	template <typename T, typename U>
		requires std::integral<T>&& std::floating_point<U>
	inline U operator()(const T& min, const U& max)
	{
		return (*this)(1.0) * (max - min) + min;
	}

	template <typename T>
		requires std::integral<T>
	inline T operator()(const T& max)
	{
		return (*this)() % (max + 1);
	}

	template <typename T>
		requires std::integral<T>
	inline T operator()(const T& min, const T& max)
	{
		return min + ((*this)() % (max - min + 1));
	}

};


class mxws_64
{

private:

	std::random_device r;

public:

	uint64_t x1, x2, w1, w2;

	typedef uint64_t result_type;

	void seed()
	{
		init();
		x1 = x2 = 1;
	}

	mxws_64(const std::seed_seq& seq)
	{
		if (seq.size() == 2)
		{
			std::vector<uint32_t> seeds(seq.size());
			seq.param(seeds.rbegin());
			w1 = (uint64_t(seeds[1]) << 32) | seeds[0];
			w2 = w1 + 1;
			x1 = x2 = 1;
		}

		else init();
	}

	mxws_64(const uint64_t& seed)
	{
		w1 = (uint64_t(r()) << 32) | r();
		w2 = seed + w1 + 1;
		x1 = x2 = 1;
	}

	void init(const uint64_t& seed)
	{
		w1 = (uint64_t(r()) << 32) | r();
		w2 = seed + w1 + 1;
		x1 = x2 = 1;
	}

	mxws_64()
	{
		init();
	}

	void init()
	{
		w1 = (uint64_t(r()) << 32) | r();
		w2 = w1 + 1;
		x1 = x2 = 1;
	}

	virtual ~mxws_64() = default;

	static uint64_t min() { return std::numeric_limits<uint64_t>::min(); }
	static uint64_t max() { return std::numeric_limits<uint64_t>::max(); }


	inline uint64_t operator()()
	{
		x1 *= w1;
		x1 = (x1 >> 32) | (x1 << 32);
		w1 += x1;

		x2 *= w2;
		x2 = (x2 >> 32) | (x2 << 32);
		w2 += x2;

		return (x1 << 32) | uint32_t(x2);
	}

	template <typename T>
		requires std::floating_point<T>
	inline T operator()(const T& f)
	{
		return ((*this)() >> 11) / 9007199254740992. * f;
	}

	template <typename T>
		requires std::floating_point<T>
	inline T operator()(const T& min, const T& max)
	{
		return (*this)(1.0) * (max - min) + min;
	}

	template <typename T, typename U>
		requires std::integral<T>&& std::floating_point<U>
	inline U operator()(const T& min, const U& max)
	{
		return (*this)(1.0) * (max - min) + min;
	}

	template <typename T>
		requires std::integral<T>
	inline T operator()(const T& max)
	{
		return (*this)() % (max + 1);
	}

	template <typename T>
		requires std::integral<T>
	inline T operator()(const T& min, const T& max)
	{
		return min + ((*this)() % (max - min + 1));
	}

};







