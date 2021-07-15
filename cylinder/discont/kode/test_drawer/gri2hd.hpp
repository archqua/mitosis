#ifndef GRI2HD_HPP
#define GRI2HD_HPP
#include <vector>

/* 2.5 dimension grid
 * that is actually just 2D grid */
template<class Value>
class Gri2hd{
	using Line = std::vector<Value>;
	std::vector<Line> grid;
public:
	Gri2hd() = default;
	Gri2hd(const Gri2hd& other) = default;
	Gri2hd(Gri2hd&& other) = default;
	Gri2hd& operator= (const Gri2hd& other) = default;
	Gri2hd& operator= (Gri2hd&& other) = default;
	Gri2hd(std::size_t I, std::size_t J){
		grid.reserve(I);
		for (std::size_t i = 0; i < I; ++i){
			grid.emplace_back(J);
		}
	}
	Line& operator[] (std::size_t i){
		return grid[i];
	}
	const Line& operator[] (std::size_t i) const {
		return grid[i];
	}
        std::size_t size() const {
          return grid.size();
        }
};

#endif
