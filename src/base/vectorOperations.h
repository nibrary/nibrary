#pragma once

#include "base/nibr.h"
#include "base/verbose.h"
#include <algorithm>
#include <vector>
#include <cmath>

namespace NIBR 
{

    template<class T1,class T2>
    void removeIdx(std::vector<T1>& inp, const std::vector<T2>& rmIdx) {
        
        if (rmIdx.size()==0) {
            return;
        } else if (inp.size()==rmIdx.size()) {
            std::vector<T1>().swap(inp);
            return;
        } 

        std::vector<T2> idx = rmIdx;

        std::sort(idx.begin(),idx.end());

        auto   beg    = inp.begin();
        std::size_t offset = 0;

        for (auto it = idx.begin(); it < idx.end(); it++) {
            std::size_t next = (it + 1 == idx.cend() ? inp.size() : *(it + 1));
            std::move(beg+*it+1, beg+next, beg+*it-offset);
            offset++;
        }
        
        inp.resize(inp.size()-idx.size());

    }

    template<class T>
    bool isUnique(std::vector<T>& inp, int ind) {

        int n = 0;

        for (int i=0; i<int(inp.size()); i++) {
            if (inp[ind] == inp[i])
                n++;
        }

        if (n==1)
            return true;
        else
            return false;

    }

    template<typename T>
    std::vector<T> linspace(T start, T end, std::size_t N) {
        std::vector<T> result;

        if (N == 0) return result;

        if (N == 1) {
            result.push_back(start);
            return result;
        }

        T step = (end - start) / static_cast<T>(N - 1);

        for (std::size_t i = 0; i < N; ++i) {
            result.push_back(start + i * step);
        }

        return result;
    }


    template<typename T>
    std::vector<T> getEvenlySeparatedSamples(const std::vector<T>& vec, std::size_t M) {
        std::vector<T> samples;
        std::size_t N = vec.size();
        
        double step = double(N-1)/double(M-1);

        for (std::size_t i = 0; i < M; ++i) {
            std::size_t index = std::size_t(std::round(i * step));
            if (index < N) {
                samples.push_back(vec[index]);
            }
        }
        
        return samples;
    }

    template<typename T>
    std::size_t findClosestIndex(const std::vector<T>& vec, const T& query) {
        if (vec.empty()) {
            disp(MSG_ERROR,"Vector is empty. Cannot find closest element.");
            return std::numeric_limits<std::size_t>::max(); // return an invalid index
        }

        std::size_t closestIndex = 0;
        T minDifference = std::numeric_limits<T>::max();

        for (std::size_t i = 0; i < vec.size(); ++i) {
            T difference = std::abs(vec[i] - query);
            if (difference < minDifference) {
                minDifference = difference;
                closestIndex = i;
            }
        }

        return closestIndex;
    }

    template<typename T>
    std::size_t findFirstGreaterIndex(const std::vector<T>& vec, const T& query) {
        if (vec.empty()) {
            disp(MSG_ERROR, "Vector is empty. Cannot find element greater than input.");
            return std::numeric_limits<std::size_t>::max(); // return an invalid index
        }

        for (std::size_t i = 0; i < vec.size(); ++i) {
            if (vec[i] >= query) {
                return i; 
            }
        }

        return std::numeric_limits<std::size_t>::max();  // return an invalid index
    }


    template<typename T>
    std::size_t findFirstSmallerIndex(const std::vector<T>& vec, const T& query) {
        if (vec.empty()) {
            disp(MSG_ERROR, "Vector is empty. Cannot find element smaller than input.");
            return std::numeric_limits<std::size_t>::max(); // return an invalid index
        }

        for (std::size_t i = 0; i < vec.size(); ++i) {
            if (vec[i] <= query) {
                return i;
            }
        }

        return std::numeric_limits<std::size_t>::max();  // return an invalid index
    }




    template <typename T>
    auto create_histogram(const std::vector<T>& data, int bin_count)
        -> std::pair<std::vector<T>, std::vector<int>> {
        if (data.empty() || bin_count <= 0) {
            return {};
        }

        const auto [min_it, max_it] = std::minmax_element(data.begin(), data.end());
        const T min_val = *min_it;
        const T max_val = *max_it;

        // Add a small epsilon if the range is zero to prevent division by zero.
        // This correctly handles the case where all data points are identical.
        double range = static_cast<double>(max_val) - static_cast<double>(min_val);
        if (range == 0.0) {
            range = 1.0; // Use a nominal range of 1 to ensure bin_width is non-zero.
        }

        const double bin_width = range / bin_count;

        std::vector<T>   bin_centers(bin_count);
        std::vector<int> bin_counts(bin_count, 0);

        // Calculate bin centers.
        for (int i = 0; i < bin_count; ++i) {
            double center_val = static_cast<double>(min_val) + (i + 0.5) * bin_width;
            // For integer types, round to the nearest value instead of truncating to improve accuracy.
            if (std::is_integral<T>::value) {
                bin_centers[i] = static_cast<T>(std::round(center_val));
            } else {
                bin_centers[i] = static_cast<T>(center_val);
            }
        }

        // Assign data points to bins.
        for (const T& value : data) {
            // Calculate bin index. The logic is robust: max_val will not exceed the bounds.
            double normalized_value = static_cast<double>(value) - static_cast<double>(min_val);
            int bin_index = static_cast<int>(normalized_value / bin_width);

            // Clamp the index to the valid range [0, bin_count - 1].
            // This handles both the max_val case (which might compute as bin_count) and
            // any floating-point inaccuracies near the upper boundary.
            bin_index = std::max(0, std::min(bin_count - 1, bin_index));

            bin_counts[bin_index]++;
        }

        return {bin_centers, bin_counts};
    }

    
    template <typename T>
    std::vector<int> create_histogram(const std::vector<T>& data, const std::vector<T>& bin_centers_const) {
        if (data.empty() || bin_centers_const.empty()) {
            return {};
        }

        // Create a sorted copy of bin centers to enable efficient searching.
        auto bin_centers = bin_centers_const;
        std::sort(bin_centers.begin(), bin_centers.end());

        std::vector<int> bin_counts(bin_centers.size(), 0);

        for (const T& value : data) {
            // Find the first bin center that is not less than the value.
            auto it = std::lower_bound(bin_centers.begin(), bin_centers.end(), value);

            int closest_bin_index = -1;

            if (it == bin_centers.begin()) {
                // Value is less than or equal to the first bin center.
                closest_bin_index = 0;
            } else if (it == bin_centers.end()) {
                // Value is greater than the last bin center.
                closest_bin_index = bin_centers.size() - 1;
            } else {
                // The value is between two centers: *(it - 1) and *it.
                // Check which one is closer.
                double dist1 = std::abs(static_cast<double>(value) - static_cast<double>(*(it - 1)));
                double dist2 = std::abs(static_cast<double>(value) - static_cast<double>(*it));

                // Note: In a tie, the lower index bin is chosen.
                if (dist1 <= dist2) {
                    closest_bin_index = std::distance(bin_centers.begin(), it - 1);
                } else {
                    closest_bin_index = std::distance(bin_centers.begin(), it);
                }
            }
            bin_counts[closest_bin_index]++;
        }

        return bin_counts;
    }

    template <typename T>
    void plot_histogram(const std::string& title, const std::vector<T>& bin_centers, const std::vector<int>& bin_counts, int max_y_axis_height = 10) {

        std::cout << "\033[0;32m" << std::flush;

        if (bin_centers.empty() || bin_counts.empty() || bin_centers.size() != bin_counts.size()) {
            std::cout << "Cannot plot histogram: Invalid input data." << std::endl;
            return;
        }

        const int max_count = *std::max_element(bin_counts.begin(), bin_counts.end());

        std::cout << std::endl << title << std::endl << std::endl;

        if (max_count == 0) {
            std::cout << "Histogram is empty (all bin counts are zero)." << std::endl;
            return;
        }

        // --- FIX: Adapt plot height for small counts ---
        // If max_count is smaller than the requested height, use max_count as the effective height.
        // This creates a 1-to-1 scale and prevents distorted Y-axis labels.
        const int effective_plot_height = std::min(max_y_axis_height, max_count);

        // --- Dynamic Width Calculation ---
        int label_width = 0;
        std::vector<std::string> labels;
        for (const T& val : bin_centers) {
            std::ostringstream oss;
            // Use a general but clean format for labels.
            if (std::abs(static_cast<double>(val)) < 1000 && std::abs(static_cast<double>(val)) > 0.1 && std::floor(val) != val) {
                oss << std::fixed << std::setprecision(1) << val;
            } else {
                oss << std::fixed << std::setprecision(0) << val;
            }
            labels.push_back(oss.str());
            if (static_cast<int>(oss.str().length()) > label_width) {
                label_width = oss.str().length();
            }
        }
        label_width = std::max(3, label_width);

        // --- Define Plotting Characters ---
        std::ostringstream bar_builder, line_builder;
        for (int i = 0; i < label_width; ++i) {
            bar_builder << "█";
            line_builder << "─";
        }
        const std::string BAR_SEGMENT = bar_builder.str();
        const std::string X_AXIS_LINE_SEGMENT = line_builder.str();
        const std::string Y_AXIS_LINE = "│";
        const std::string X_AXIS_TICK_DOWN = "┴";
        const std::string AXIS_CORNER = "└";
        const std::string EMPTY_SLOT(label_width, ' ');
        
        const double scale = static_cast<double>(max_count) / effective_plot_height;
        const int y_label_interval = std::max(1, effective_plot_height / 4);

        // Plot Y-axis and bars using the effective_plot_height
        for (int y_level = effective_plot_height; y_level > 0; --y_level) {
            if (y_level == effective_plot_height || (y_level % y_label_interval == 0)) {
                std::cout << std::setw(9) << std::fixed << std::setprecision(0) << std::ceil(y_level * scale) << " " << Y_AXIS_LINE;
            } else {
                std::cout << "          " << Y_AXIS_LINE;
            }

            for (size_t i = 0; i < bin_counts.size(); ++i) {
                std::cout << " ";
                // The check to draw a bar segment is now more intuitive
                if (static_cast<double>(bin_counts[i]) >= y_level * scale) {
                    std::cout << BAR_SEGMENT;
                } else {
                    std::cout << EMPTY_SLOT;
                }
            }
            std::cout << std::endl;
        }

        // Plot X-axis line
        std::cout << "        0 " << AXIS_CORNER << "─";
        for (size_t i = 0; i < bin_counts.size(); ++i) {
            std::cout << X_AXIS_LINE_SEGMENT << (i < bin_counts.size() - 1 ? X_AXIS_TICK_DOWN : "─");
        }
        std::cout << std::endl;

        // Plot X-axis labels with dynamic padding
        std::cout << "           "; // Align with content
        for (size_t i = 0; i < labels.size(); ++i) {
            int padding_total = label_width - labels[i].length();
            int padding_left  = padding_total / 2;
            int padding_right = padding_total - padding_left;
            std::cout << " " << std::string(padding_left, ' ') << labels[i] << std::string(padding_right, ' ');
        }
        std::cout << std::endl << std::endl;

        std::cout << "\033[0m" << std::flush;
    }
    
}
