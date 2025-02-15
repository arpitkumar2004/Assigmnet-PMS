#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <string>
#include <functional>
#include <algorithm>

class SeparationColumn {
public:
    SeparationColumn(std::map<std::string, double> input_feed, std::map<std::string, double> output_vapor,
                     std::map<std::string, double> output_liquid, std::map<std::string, double> volatility_factors)
        : feed(input_feed), vapor(output_vapor), liquid(output_liquid), volatilities(volatility_factors) {
        relative_volatility = compute_relative_volatility(volatilities["iC5"], volatilities["nC5"]);
    }

    void execute() {
        std::cout << "Calculated relative volatility: " << relative_volatility << std::endl;
        double min_stages = compute_min_stages();
        std::cout << "Required minimum stages: " << min_stages << std::endl;

        std::cout << "Distribution of non-key components:" << std::endl;
        auto nonkey_results = compute_nonkey_distribution();
        for (const auto& comp : nonkey_results) {
            std::cout << comp.first << " -> Vapor: " << comp.second.first << ", Liquid: " << comp.second.second << std::endl;
        }

        double min_reflux_ratio = compute_min_reflux_ratio();
        std::cout << "Calculated minimum reflux ratio: " << min_reflux_ratio << std::endl;
        double actual_stages = estimate_actual_stages(min_reflux_ratio, min_stages);
        std::cout << "Estimated actual stages: " << actual_stages << std::endl;
    }

private:
    std::map<std::string, double> feed;
    std::map<std::string, double> vapor;
    std::map<std::string, double> liquid;
    std::map<std::string, double> volatilities;
    double relative_volatility;

    double compute_relative_volatility(double key_high, double key_low) {
        return key_high / key_low;
    }

    double compute_min_stages() {
        double vapor_HK = vapor["iC5"] / sum_composition(vapor);
        double liquid_HK = liquid["iC5"] / sum_composition(liquid);
        double vapor_LK = vapor["nC5"] / sum_composition(vapor);
        double liquid_LK = liquid["nC5"] / sum_composition(liquid);
        return log((vapor_HK / vapor_LK) * (liquid_LK / liquid_HK)) / log(relative_volatility);
    }

    std::map<std::string, std::pair<double, double>> compute_nonkey_distribution() {
        std::map<std::string, std::pair<double, double>> result;
        for (const auto& comp : feed) {
            if (comp.first != "iC5" && comp.first != "nC5") {
                result[comp.first] = {vapor[comp.first], liquid[comp.first]};
            }
        }
        return result;
    }

    double sum_composition(const std::map<std::string, double>& composition) {
        double sum = 0.0;
        for (const auto& pair : composition) {
            sum += pair.second;
        }
        return sum;
    }

    double find_theta(double q = 1.0) {
        double total_feed = sum_composition(feed);
        std::vector<double> alpha_values;
        std::vector<double> z_values;
        for (const auto& comp : volatilities) {
            alpha_values.push_back(comp.second / volatilities["nC5"]);
            z_values.push_back(feed[comp.first] / total_feed);
        }
        
        double theta_low = 0.0, theta_high = 2.0, theta_mid;
        double tolerance = 1e-6;
        while (theta_high - theta_low > tolerance) {
            theta_mid = (theta_low + theta_high) / 2.0;
            double sum = 0.0;
            for (size_t i = 0; i < alpha_values.size(); ++i) {
                sum += alpha_values[i] * z_values[i] / (alpha_values[i] - theta_mid);
            }
            sum -= (1 - q);
            if (sum > 0) {
                theta_high = theta_mid;
            } else {
                theta_low = theta_mid;
            }
        }
        return theta_mid;
    }

    double compute_min_reflux_ratio() {
        double theta = find_theta();
        double total_vapor = sum_composition(vapor);
        double sum = 0.0;
        for (const auto& comp : volatilities) {
            double vapor_frac = vapor[comp.first] / total_vapor;
            double alpha_i = comp.second / volatilities["nC5"];
            sum += alpha_i * vapor_frac / (alpha_i - theta);
        }
        return sum - 1;
    }

    double estimate_actual_stages(double min_reflux_ratio, double min_stages) {
        double reflux_ratio = 1.3 * min_reflux_ratio;
        double x = (reflux_ratio - min_reflux_ratio) / (reflux_ratio + 1);
        double y = 1 - exp(((1 + 54.4 * x) * (x - 1)) / ((11 + 117.2 * x) * sqrt(x)));
        return (min_stages + 1) / (1 - y);
    }
};

int main() {
    std::map<std::string, double> input_feed = { {"iC4", 12}, {"nC4", 448}, {"iC5", 36}, {"nC5", 15}, {"C6", 23}, {"C7", 39.1}, {"C8", 272.2}, {"C9", 31.0} };
    std::map<std::string, double> output_vapor = { {"iC4", 12}, {"nC4", 442}, {"iC5", 13}, {"nC5", 1} };
    std::map<std::string, double> output_liquid = { {"nC4", 6}, {"iC5", 23}, {"nC5", 14}, {"C6", 23}, {"C7", 39.1}, {"C8", 272.2}, {"C9", 31.0} };
    std::map<std::string, double> volatility_factors = { {"iC4", 3.0}, {"nC4", 2.5}, {"iC5", 1.2}, {"nC5", 1.0}, {"C6", 0.6}, {"C7", 0.25}, {"C8", 0.12}, {"C9", 0.1} };

    SeparationColumn column(input_feed, output_vapor, output_liquid, volatility_factors);
    column.execute();
    return 0;
}


