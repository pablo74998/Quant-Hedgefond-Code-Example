// quant_fund.cpp
// Beispielimplementierung eines quantitativen Hedgefonds in C++
// Enthält Data Loading, Feature Engineering, Label Generation,
// Strategy Model (Placeholder), Risk Management, Slippage & TCA,
// Portfolio Optimizer, Multi-Asset Support und Live Execution Stub.

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <random>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <chrono>
#include <sstream>
#include <iomanip>

// === Typedefs ===
using Date = std::string;                        // YYYY-MM-DD
using PriceSeries = std::map<Date, double>;
using MultiPrice = std::map<std::string, PriceSeries>;
using Weights = std::map<std::string, double>;

// === Data Loader ===
class DataLoader {
private:
    int days_in_month(int year, int month) {
        if (month == 2) {
            return 28; // Ignoriert Schaltjahre für Vereinfachung
        } else if (month == 4 || month == 6 || month == 9 || month == 11) {
            return 30;
        } else {
            return 31;
        }
    }

public:
    DataLoader(const std::vector<std::string>& tickers, const Date& start, const Date& end)
        : tickers_(tickers), start_(start), end_(end) {}

    MultiPrice fetch_prices() {
        MultiPrice data;
        std::mt19937_64 rng(42);
        std::normal_distribution<double> dist(0, 1);
        int days = 200;

        // Generiere Datumsreihe
        std::vector<Date> dates;
        int year = 2020;
        int month = 1;
        int day = 1;
        for (int i = 0; i < days; ++i) {
            std::ostringstream oss;
            oss << year << "-" 
                << std::setw(2) << std::setfill('0') << month << "-"
                << std::setw(2) << std::setfill('0') << day;
            dates.push_back(oss.str());

            // Tagesinkrement
            day++;
            if (day > days_in_month(year, month)) {
                day = 1;
                month++;
                if (month > 12) {
                    month = 1;
                    year++;
                }
            }
        }

        for (const auto& t : tickers_) {
            double price = 100.0;
            PriceSeries series;
            for (const auto& d : dates) {
                price += dist(rng);
                series[d] = price;
            }
            data[t] = series;
        }
        return data;
    }

private:
    std::vector<std::string> tickers_;
    Date start_, end_;
};

// === Feature Engineering ===
struct Features {
    std::map<std::string, double> momentum;
    std::map<std::string, double> volatility;
};

class FeatureEngineer {
public:
    FeatureEngineer(int mom_window = 5, int vol_window = 20)
        : mom_w(mom_window), vol_w(vol_window) {}

    std::map<Date, Features> compute(const MultiPrice& prices) {
        std::map<Date, Features> feats;
        for (const auto& [ticker, series] : prices) {
            std::vector<double> vals;
            std::vector<Date> dates;
            for (const auto& [date, price] : series) {
                dates.push_back(date);
                vals.push_back(price);
            }
            int n = vals.size();
            for (int i = vol_w; i < n; ++i) {
                // Momentum berechnen
                double mom = (vals[i] - vals[i - mom_w]) / vals[i - mom_w];
                // Volatilität berechnen
                double sum = 0.0, sq = 0.0;
                for (int j = i - vol_w + 1; j <= i; ++j) {
                    sum += vals[j];
                    sq += vals[j] * vals[j];
                }
                double mean = sum / vol_w;
                double var = (sq / vol_w) - (mean * mean);
                var = std::max(var, 0.0); // Vermeidet negative Varianz
                double vol = std::sqrt(var);
                feats[dates[i]].momentum[ticker] = mom;
                feats[dates[i]].volatility[ticker] = vol;
            }
        }
        return feats;
    }

private:
    int mom_w, vol_w;
};

// === Label Generation ===
class LabelGenerator {
public:
    LabelGenerator(double thresh = 0.005) : threshold(thresh) {}

    std::map<Date, Weights> generate(const MultiPrice& prices) {
        std::map<Date, Weights> labels;
        for (const auto& [ticker, series] : prices) {
            std::vector<Date> dates;
            std::vector<double> vals;
            for (const auto& [date, price] : series) {
                dates.push_back(date);
                vals.push_back(price);
            }
            for (size_t i = 0; i < vals.size() - 1; ++i) {
                double ret = (vals[i + 1] - vals[i]) / vals[i];
                labels[dates[i]][ticker] = (ret > threshold) ? 1.0 : 0.0;
            }
        }
        return labels;
    }

private:
    double threshold;
};

// === Risk Management (angepasst) ===
class RiskManager {
public:
    RiskManager(double max_dd = 0.1, double max_pos = 0.2)
        : max_drawdown(max_dd), max_position(max_pos), peak(0) {}

    Weights size_positions(const Weights& signal,
                           const std::map<std::string, std::vector<double>>& returns) {
        Weights weighted;
        std::map<std::string, double> inv_vol;

        for (const auto& [ticker, rets] : returns) {
            if (rets.empty()) {
                inv_vol[ticker] = 0.0;
                continue;
            }
            double sum = std::accumulate(rets.begin(), rets.end(), 0.0);
            double mean = sum / rets.size();
            double sq_sum = std::inner_product(rets.begin(), rets.end(), rets.begin(), 0.0);
            double var = (sq_sum / rets.size()) - (mean * mean);
            var = std::max(var, 0.0);
            double vol = std::sqrt(var);
            inv_vol[ticker] = (vol > 1e-6) ? 1.0 / vol : 0.0;
        }

        double total = 0.0;
        for (const auto& [ticker, s] : signal) {
            weighted[ticker] = s * inv_vol[ticker];
            total += std::abs(weighted[ticker]);
        }

        if (total > 1e-6) {
            for (auto& [ticker, w] : weighted) {
                w /= total;
                // Manuelle Clamp-Implementierung für C++11-Kompatibilität
                w = std::max(-max_position, std::min(w, max_position));
            }
        } else {
            for (auto& [ticker, w] : weighted) {
                w = 0.0;
            }
        }
        return weighted;
    }

    // ... (andere Methoden unverändert)

private:
    double max_drawdown, max_position, peak;
};

// === Backtest-Loop mit korrekten Renditen ===
int main() {
    std::vector<std::string> tickers = {"AAPL", "MSFT", "GOOG", "AMZN"};
    DataLoader loader(tickers, "2020-01-01", "2020-10-01");
    MultiPrice prices = loader.fetch_prices();

    FeatureEngineer fe;
    auto features = fe.compute(prices);

    StrategyModel model;
    RiskManager risk(0.15, 0.25);
    CostModel cost_model;
    PortfolioOptimizer optimizer;
    LiveExecutor live("API_KEY", "SECRET", "https://api.example.com");

    double equity = 1e6;
    std::map<std::string, std::vector<double>> return_history;

    for (const auto& [date, feat] : features) {
        Weights signal = model.predict(feat);
        Weights positions = risk.size_positions(signal, return_history);
        positions = optimizer.optimize(positions);

        // Renditen basierend auf nächstem Preis berechnen
        std::map<std::string, double> daily_returns;
        for (const auto& ticker : tickers) {
            auto& series = prices[ticker];
            auto current = series.find(date);
            if (current == series.end()) continue;
            auto next = std::next(current);
            if (next == series.end()) continue;

            double ret = (next->second - current->second) / current->second;
            daily_returns[ticker] = ret;
            return_history[ticker].push_back(ret);
        }

        // Portfoliorendite und Kosten anwenden
        double portfolio_ret = 0.0;
        for (const auto& [ticker, weight] : positions) {
            portfolio_ret += weight * daily_returns[ticker];
            equity -= cost_model.apply(ticker, weight, equity);
        }
        equity *= (1 + portfolio_ret);

        // Risikomanagement-Update
        risk.update_peak(equity);
        if (!risk.check_drawdown(equity)) {
            std::cout << "Drawdown-Limit überschritten. Handel gestoppt.\n";
            break;
        }
    }

    std::cout << "Endkapital: $" << equity << std::endl;
    return 0;
}