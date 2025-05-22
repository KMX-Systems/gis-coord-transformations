/// Copyright (c) 2025 - present KMX Systems. All rights reserved.
/// @file gis/uniform_grid_mapper.hpp
/// @brief Defines a 2D grid.
///
/// @note This file assumes a C++20 compliant compiler.
/// @warning Does not use PCH if PCH is not defined.
#pragma once
#ifndef PCH
    #include <kmx/gis/coordinate/base.hpp>
#endif

/// @namespace kmx::gis
/// @brief Namespace for Geographic Information System (GIS) related functionalities.
///
/// This namespace contains structures, constants, and algorithms commonly used
/// in geodetic and cartographic computations.
namespace kmx::gis
{
    /// @brief Maps points to cell indices in a uniform 2D grid.
    ///
    /// This class calculates the 0-based row and column index for a given 2D point
    /// within a defined rectangular grid. The grid is characterized by its minimum
    /// and maximum corners and the number of rows and columns.
    /// The cell indices are returned as a `kmx::gis::coordinate::xy` object where `xy.x` is the column
    /// and `xy.y` is the row.
    ///
    /// @tparam PointType_ The type of 2D point to work with (e.g., `kmx::gis::coordinate::xy<double>`).
    ///                    It must provide public members `x` and `y`, and a `value_type` alias
    ///                    for its coordinate type (which should typically be floating-point).
    ///                    Defaults to `kmx::gis::coordinate::xy<double>`.
    /// @tparam IndexComponentType_ The unsigned integral type for the components of the cell identifier
    ///                             (e.g., `std::uint8_t`, `unsigned int`). Defaults to `std::uint8_t`.
    /// @tparam IntermediateFPType_ The floating-point type used for internal normalization and scaling calculations
    ///                             (e.g., `float`, `double`). Defaults to `double`.
    template <typename PointType_ = coordinate::xy<double>, typename IndexComponentType_ = std::uint8_t,
              typename IntermediateFPType_ = double>
    class uniform_grid_mapper
    {
        /// @brief Static assertion to ensure PointType_::value_type is arithmetic (typically float for spatial points).
        static_assert(std::is_arithmetic_v<typename PointType_::value_type>, "PointType_::value_type must be an arithmetic type.");
        /// @brief Static assertion to ensure IndexComponentType_ is an unsigned integral type.
        static_assert(std::is_integral_v<IndexComponentType_> && std::is_unsigned_v<IndexComponentType_>,
                      "IndexComponentType_ for cell components must be an unsigned integral type.");
        /// @brief Static assertion to ensure IntermediateFPType_ is a floating-point type.
        static_assert(std::is_floating_point_v<IntermediateFPType_>, "IntermediateFPType_ must be a floating-point type.");

    public:
        using point_type = PointType_;                           ///< Alias for the point type used by the mapper.
        using coordinate_type = typename point_type::value_type; ///< Alias for the PointType_'s coordinate underlying type.
        using index_component_type = IndexComponentType_;        ///< Alias for the type of each component of a cell index.
        using intermediate_fp_type = IntermediateFPType_;        ///< Alias for the internal calculation FP type.

        /// @brief Alias for the grid cell identifier, using `kmx::gis::coordinate::xy<IndexComponentType_>`.
        /// Its `.x` member represents the column index, and its `.y` member represents the row index.
        using cell_type = kmx::gis::coordinate::xy<index_component_type>;

        static constexpr auto zero_index = static_cast<index_component_type>(0);
        static constexpr auto zero_value = static_cast<coordinate_type>(0);

        /// @brief Constructs a uniform_grid_mapper.
        /// @param min_corner The minimum corner (e.g., bottom-left) of the grid's spatial extent.
        /// @param max_corner The maximum corner (e.g., top-right) of the grid's spatial extent.
        /// @param row_count The number of rows in the grid. Must be positive.
        /// @param column_count The number of columns in the grid. Must be positive.
        /// @throws std::invalid_argument If row_count, column_count are not positive, or if the extent has non-positive width/height.
        constexpr uniform_grid_mapper(const point_type& min_corner, const point_type& max_corner, const index_component_type row_count,
                                      const index_component_type column_count) noexcept(false):
            min_corner_(min_corner),
            max_corner_(max_corner),
            row_count_(row_count),
            column_count_(column_count),
            total_width_(max_corner_.x - min_corner_.x),
            total_height_(max_corner_.y - min_corner_.y)
        {
            if (!std::is_constant_evaluated())
            {
                if ((row_count_ == zero_index) || (column_count_ == zero_index))
                    throw std::invalid_argument("Number of rows and columns must be positive.");

                if ((total_width_ <= zero_value) || (total_height_ <= zero_value))
                    throw std::invalid_argument("Spatial extent must have positive width and height.");
            }
        }

        /// @brief Maps a given point to its corresponding grid cell indices using array-like access.
        /// @param point The 2D point to map.
        /// @return An std::optional containing the cell_type (xy.x = column, xy.y = row)
        ///         if the point is within the grid extent, otherwise an empty std::optional.
        [[nodiscard]] constexpr std::optional<cell_type> operator[](const point_type& point) const noexcept
        {
            /// Epsilon for boundary checks, adaptive to coordinate_type's precision.
            constexpr coordinate_type boundary_epsilon = []
            {
                if constexpr (std::is_same_v<coordinate_type, float>)
                    return 1e-6f;
                else if constexpr (std::is_same_v<coordinate_type, double>)
                    return 1e-9;
                else if constexpr (std::is_same_v<coordinate_type, long double>)
                    return 1e-12L;
                else
                    return static_cast<coordinate_type>(0);
            }();

            // Check if point is outside the grid extent (with epsilon tolerance)
            if ((point.x < (min_corner_.x - boundary_epsilon)) || (point.x > (max_corner_.x + boundary_epsilon)) ||
                (point.y < (min_corner_.y - boundary_epsilon)) || (point.y > (max_corner_.y + boundary_epsilon)))
                return {};

            // Clamp coordinates to the strict grid boundaries
            const coordinate_type clamped_x = std::max(min_corner_.x, std::min(point.x, max_corner_.x));
            const coordinate_type clamped_y = std::max(min_corner_.y, std::min(point.y, max_corner_.y));

            // Robustness checks for dimensions
            if ((total_width_ <= static_cast<coordinate_type>(0)) || (total_height_ <= static_cast<coordinate_type>(0)))
                return {};
            if ((column_count_ == static_cast<index_component_type>(0)) || (row_count_ == static_cast<index_component_type>(0)))
                return {};

            // Normalize coordinates to [0, 1] using IntermediateFPType_
            const intermediate_fp_type norm_x =
                static_cast<intermediate_fp_type>(clamped_x - min_corner_.x) / static_cast<intermediate_fp_type>(total_width_);
            const intermediate_fp_type norm_y =
                static_cast<intermediate_fp_type>(clamped_y - min_corner_.y) / static_cast<intermediate_fp_type>(total_height_);

            // Scale normalized coordinates
            const intermediate_fp_type scaled_col = norm_x * static_cast<intermediate_fp_type>(column_count_);
            const intermediate_fp_type scaled_row = norm_y * static_cast<intermediate_fp_type>(row_count_);

            // Floor to get raw indices
            intermediate_fp_type floored_col = std::floor(scaled_col);
            intermediate_fp_type floored_row = std::floor(scaled_row);

            // Zero constant
            constexpr auto zero = static_cast<intermediate_fp_type>(0);

            // Safeguard against negative results from floor
            if (floored_col < zero)
                floored_col = zero;
            if (floored_row < zero)
                floored_row = zero;

            // Clamp floored values to the maximum representable by IndexComponentType_ before casting
            constexpr intermediate_fp_type max_index_fp_val =
                static_cast<intermediate_fp_type>(std::numeric_limits<index_component_type>::max());
            if (floored_col > max_index_fp_val)
                floored_col = max_index_fp_val;
            if (floored_row > max_index_fp_val)
                floored_row = max_index_fp_val;

            // Cast to final IndexComponentType_
            index_component_type col_idx = static_cast<index_component_type>(floored_col);
            index_component_type row_idx = static_cast<index_component_type>(floored_row);

            // Ensure indices are within the valid grid cell range [0, N-1]
            col_idx = std::min(col_idx, static_cast<index_component_type>(column_count_ - 1));
            row_idx = std::min(row_idx, static_cast<index_component_type>(row_count_ - 1));

            // Construct cell_type (which is xy<IndexComponentType_>) with (x=col_idx, y=row_idx)
            return cell_type {col_idx, row_idx};
        }

        /// @brief Gets the dimensions (row count, column count) of the grid.
        /// @return A std::pair of IndexComponentType_ for (row_count, column_count).
        [[nodiscard]] constexpr std::pair<index_component_type, index_component_type> size() const noexcept
        {
            return {row_count_, column_count_};
        }

    private:
        /// @brief The minimum corner (e.g., bottom-left) of the grid's spatial extent.
        point_type min_corner_;
        /// @brief The maximum corner (e.g., top-right) of the grid's spatial extent.
        point_type max_corner_;
        /// @brief The number of rows in the grid.
        index_component_type row_count_;
        /// @brief The number of columns in the grid.
        index_component_type column_count_;

        /// @brief The total width of the grid's spatial extent.
        coordinate_type total_width_;
        /// @brief The total height of the grid's spatial extent.
        coordinate_type total_height_;
    };

} // namespace kmx::gis
