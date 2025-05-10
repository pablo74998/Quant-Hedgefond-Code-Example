// ... (vorherige Includes bleiben unverändert)

// === Erweiterter Risk Manager mit Risikoanalyse ===
class RiskManager {
private:
    double max_drawdown, max_position, peak;
    double max_sector_exposure;       // Maximaler Sektoranteil
    double var_threshold;             // VaR-Schwellenwert
    double stress_loss_threshold;     // Stress-Test-Limit
    std::vector<double> portfolio_returns; // Portfolio-Renditenverlauf
    std::map<std::string, std::string> ticker_to_sector; // Sektor-Zuordnung

public:
    RiskManager(double max_dd = 0.1, double max_pos = 0.2, 
                double max_sector = 0.5, double var_thresh = -0.05,
                double stress_loss = -0.1)
        : max_drawdown(max_dd), max_position(max_pos), 
          max_sector_exposure(max_sector), var_threshold(var_thresh),
          stress_loss_threshold(stress_loss) {
        // Beispielhafte Sektorzuordnung
        ticker_to_sector = {
            {"AAPL", "Tech"}, {"MSFT", "Tech"}, 
            {"GOOG", "Tech"}, {"AMZN", "Tech"}
        };
    }

    // Position Sizing mit Volatilitätsgewichtung (bestehend)
    Weights size_positions(const Weights& signal,
                          const std::map<std::string, std::vector<double>>& returns) {
        // ... (unverändert vom Originalcode)
    }

    // Risiko-Checks für Positionsgrößen
    Weights apply_risk_checks(const Weights& positions,
                             const std::map<std::string, std::vector<double>>& asset_returns) {
        Weights adjusted = positions;

        // 1. Sektor-Exposure prüfen
        std::map<std::string, double> sector_exposure;
        for (const auto& [ticker, weight] : adjusted) {
            sector_exposure[ticker_to_sector[ticker]] += std::abs(weight);
        }
        for (auto& [sector, exposure] : sector_exposure) {
            if (exposure > max_sector_exposure) {
                // Reduziere Positionen im Sektor linear
                for (auto& [ticker, weight] : adjusted) {
                    if (ticker_to_sector[ticker] == sector) {
                        weight *= 0.5; // Beispielhafte Reduktion
                    }
                }
            }
        }

        // 2. VaR-Check (95% Konfidenz)
        double var = calculate_VaR(0.95);
        if (var < var_threshold) {
            // Globale Positionsreduktion bei VaR-Verletzung
            for (auto& [ticker, weight] : adjusted) {
                weight *= 0.7;
            }
        }

        // 3. Stress-Test: Simuliere simultanen Crash aller Assets
        std::map<std::string, double> crash_scenario;
        for (const auto& [ticker, _] : adjusted) {
            crash_scenario[ticker] = -0.25; // 25% Verlust annehmen
        }
        double stress_loss = simulate_stress_test(adjusted, crash_scenario);
        if (stress_loss < stress_loss_threshold) {
            for (auto& [ticker, weight] : adjusted) {
                weight *= 0.6; // Aggressive Reduktion
            }
        }

        // 4. Korrelations-Check (Reduziere hochkorrelierte Positionen)
        auto corr_matrix = compute_correlations(asset_returns);
        for (const auto& [ticker1, row] : corr_matrix) {
            for (const auto& [ticker2, corr] : row) {
                if (ticker1 != ticker2 && corr > 0.7) {
                    adjusted[ticker1] *= 0.8;
                    adjusted[ticker2] *= 0.8;
                }
            }
        }

        return adjusted;
    }

    // Hilfsmethoden für Risikoanalysen
    void track_portfolio_return(double ret) { 
        portfolio_returns.push_back(ret); 
    }

    double calculate_VaR(double confidence = 0.95) {
        if (portfolio_returns.empty()) return 0.0;
        auto sorted = portfolio_returns;
        std::sort(sorted.begin(), sorted.end());
        int idx = static_cast<int>((1 - confidence) * sorted.size());
        return sorted[std::max(idx, 0)];
    }

    double calculate_ES(double confidence = 0.95) {
        double var = calculate_VaR(confidence);
        double sum = 0.0;
        int count = 0;
        for (double ret : portfolio_returns) {
            if (ret <= var) {
                sum += ret;
                count++;
            }
        }
        return (count > 0) ? sum / count : 0.0;
    }

    double simulate_stress_test(const Weights& positions,
                               const std::map<std::string, double>& scenario) {
        double loss = 0.0;
        for (const auto& [ticker, shock] : scenario) {
            if (positions.count(ticker)) {
                loss += positions.at(ticker) * shock;
            }
        }
        return loss;
    }

    std::map<std::string, std::map<std::string, double>> compute_correlations(
        const std::map<std::string, std::vector<double>>& asset_returns) {
        std::map<std::string, std::map<std::string, double>> corr_matrix;
        for (const auto& [ticker1, rets1] : asset_returns) {
            for (const auto& [ticker2, rets2] : asset_returns) {
                if (ticker1 >= ticker2) continue; // Symmetrie ausnutzen
                
                size_t n = std::min(rets1.size(), rets2.size());
                if (n < 2) {
                    corr_matrix[ticker1][ticker2] = 0.0;
                    continue;
                }

                double sum1 = 0.0, sum2 = 0.0, sum12 = 0.0;
                double sq1 = 0.0, sq2 = 0.0;
                for (size_t i = 0; i < n; ++i) {
                    sum1 += rets1[i];
                    sum2 += rets2[i];
                    sum12 += rets1[i] * rets2[i];
                    sq1 += rets1[i] * rets1[i];
                    sq2 += rets2[i] * rets2[i];
                }

                double cov = (sum12 - sum1*sum2/n) / n;
                double var1 = (sq1 - sum1*sum1/n) / n;
                double var2 = (sq2 - sum2*sum2/n) / n;
                
                double corr = (var1 > 0 && var2 > 0) ? 
                    cov / (sqrt(var1)*sqrt(var2)) : 0.0;
                
                corr_matrix[ticker1][ticker2] = corr;
                corr_matrix[ticker2][ticker1] = corr;
            }
        }
        return corr_matrix;
    }

    // Drawdown-Logik (bestehend)
    void update_peak(double equity) { 
        if (equity > peak) peak = equity; 
    }
    
    bool check_drawdown(double equity) { 
        return (peak - equity)/peak <= max_drawdown; 
    }
};

// ... (Rest des bestehenden Codes bleibt unverändert bis zur main-Funktion)

int main() {
    // Initialisierung wie zuvor
    std::vector<std::string> tickers = {"AAPL", "MSFT", "GOOG", "AMZN"};
    DataLoader loader(tickers, "2020-01-01", "2020-10-01");
    MultiPrice prices = loader.fetch_prices();

    FeatureEngineer fe;
    auto features = fe.compute(prices);

    LabelGenerator label_gen;
    auto labels = label_gen.generate(prices);

    // Risk Manager mit erweiterten Parametern
    RiskManager risk(
        0.15,    // max_drawdown
        0.25,    // max_position
        0.4,     // max_sector_exposure
        -0.03,   // var_threshold (3% Tages-VaR)
        -0.08    // stress_loss_threshold (8% Stress-Verlust)
    );

    // ... (Rest der Initialisierung)

    double equity = 1e6;
    std::map<std::string, std::vector<double>> return_history;

    for (const auto& [date, feat] : features) {
        // ... (Signalgenerierung wie zuvor)
        
        Weights positions = risk.size_positions(signal, return_history);
        positions = risk.apply_risk_checks(positions, return_history); // Neu!
        positions = optimizer.optimize(positions);

        // ... (Portfoliobewertung und Renditeberechnung)

        // Risikometriken aktualisieren
        risk.track_portfolio_return(portfolio_ret);
        
        // Risikostatus ausgeben
        if (i % 20 == 0) { // Jeden 20. Tag ausgeben
            std::cout << "Risikoreport " << date << ":\n";
            std::cout << "  VaR 95%: " << risk.calculate_VaR()*100 << "%\n";
            std::cout << "  ES 95%: " << risk.calculate_ES()*100 << "%\n";
        }

        // ... (Drawdown-Check wie zuvor)
    }

    std::cout << "Endkapital: $" << equity << std::endl;
    return 0;
}