/* Mechanical Refrigeration Cycle */
#include <stdio.h>
#include <math.h>
#define PI 3.14159265358979323846 // π = 3.14159265358979323846
#define GAMMA 1.4 // Heat capacity ratio for ideal gas 
#define THERMAL_CONDUCTIVITY_WATER 0.606 // Thermal conductivity of water in W/(m·K) at 25°C
#define DEW_POINT_A 17.27 // Empirical constant for dew point calculation
#define DEW_POINT_B 237.7 // Empirical constant for dew point calculation

// Heat Pump Modes
#define COOLING_MODE 1
#define HEATING_MODE 2

float calculateWorkDone(float h_out, float h_in, float m);
float calculateHeatRejected(float h_out, float h_in, float m, float latentHeat);
float calculateHeatAbsorbed(float h_out, float h_in, float m, float latentHeat);
float calculateCOP(float Qe, float Wc_actual);
float calculateActualCompressorWork(float h2s, float h1, float efficiency, float m);
float calculateActualTurbineWork(float h2, float h1s, float efficiency, float m);
float adjustEnthalpyForPressureDrop(float h, float deltaP);
float getSaturatedLiquidEnthalpy(float T, float P);
float getSaturatedVaporEnthalpy(float T, float P);
float getLatentHeat(float T, float P);
float calculatePressureDrop(float L, float D, float rho, float v, float f);
float calculateExpansionPressureDrop(float P1, float T1, float T2);
float getSpecificHeatCapacity(float T, float P, int isConstantPressure);
float calculateIsentropicEnthalpy(float h1, float P1, float P2);
float calculateIdealWork(float h1, float h2, float m);
float calculateSecondLawEfficiency(float W_actual, float W_ideal);
float calculateHeatTransferCoefficient(float Re, float Pr, float D, float foulingFactor);
float calculateHeatTransferRate(float h, float area, float dT);
float calculateHumidityRatio(float P_v, float P);
float calculateDewPoint(float T, float RH);
float calculateWetBulbTemperature(float T, float T_dp);

int main() {
    float L = 10.0;  // Length of pipe in evaporator and condenser (m)
    float D = 0.05;  // Diameter of pipe in evaporator and condenser (m)
    float rho = 1.0; // Density of the fluid (kg/m³) - assume a value
    float v = 2.0;   // Velocity of the fluid (m/s) - assume a value
    float Pr = 7.0;  // Prandtl number - assume a value for the refrigerant
    float T_evap = 5.0;  // Saturation temperature in the evaporator (°C)
    float T_cond = 35.0; // Saturation temperature in the condenser (°C)
    float P_evap = 1.2;  // Saturation pressure in the evaporator (bar)
    float P_cond = 4.5;  // Saturation pressure in the condenser (bar)
    float deltaP_e = calculatePressureDrop(L, D, rho, v, 0.02); // Friction factor (0.02 is an example)
    float deltaP_c = calculatePressureDrop(L, D, rho, v, 0.02); // Friction factor (0.02 is an example)
    float deltaP_v = calculateExpansionPressureDrop(P_evap, T_evap, T_cond);
    float eta_c = 0.85; // Compressor efficiency
    float eta_t = 0.90; // Turbine efficiency
    float eta_e = 0.90; // Expansion valve efficiency
    float foulingFactor_evap = 0.05; // Fouling factor for evaporator
    float foulingFactor_cond = 0.05; // Fouling factor for condenser

    // Mass flow rate (kg/s)
    float m = 1.0;

    // Heat Pump Mode Selection
    int mode = COOLING_MODE; // Set mode to COOLING_MODE or HEATING_MODE

    float h1_liquid = getSaturatedLiquidEnthalpy(T_evap, P_evap); // Enthalpy of saturated liquid in evaporator
    float h1_vapor = getSaturatedVaporEnthalpy(T_evap, P_evap);  // Enthalpy of saturated vapor in evaporator
    float h2_vapor = getSaturatedVaporEnthalpy(T_cond, P_cond);  // Enthalpy of saturated vapor at compressor outlet
    float h3_liquid = getSaturatedLiquidEnthalpy(T_cond, P_cond); // Enthalpy of saturated liquid at condenser outlet
    float h4_vapor = getSaturatedVaporEnthalpy(T_evap, P_evap);  // Enthalpy of saturated vapor at expansion valve outlet

    // Get latent heat values
    float latentHeat_evap = getLatentHeat(T_evap, P_evap);
    float latentHeat_cond = getLatentHeat(T_cond, P_cond);

    // Isentropic enthalpies
    float h2s = calculateIsentropicEnthalpy(h1_vapor, P_evap, P_cond);
    float h1s = calculateIsentropicEnthalpy(h2_vapor, P_cond, P_evap);

    float h1_adjusted = adjustEnthalpyForPressureDrop(h1_vapor, deltaP_e);
    float h2_adjusted = adjustEnthalpyForPressureDrop(h2_vapor, deltaP_c);
    float h3_adjusted = adjustEnthalpyForPressureDrop(h3_liquid, deltaP_c);
    float h4_adjusted = adjustEnthalpyForPressureDrop(h4_vapor, deltaP_v);

    // Actual work done by the compressor
    float Wc_actual = calculateActualCompressorWork(h2_vapor, h1_vapor, eta_c, m);

    // Actual work done by the turbine 
    float Wt_actual = calculateActualTurbineWork(h1s, h2_vapor, eta_t, m);

    // Ideal work (isentropic work) for comparison
    float Wc_ideal = calculateIdealWork(h2s, h1_vapor, m);
    float Wt_ideal = calculateIdealWork(h2_vapor, h1s, m);

    // Heat rejected or absorbed depending on mode
    float Qc, Qe;

    if (mode == COOLING_MODE) {
        Qc = calculateHeatRejected(h2_adjusted, h3_adjusted, m, latentHeat_cond);
        Qe = calculateHeatAbsorbed(h1_adjusted, h4_adjusted, m, latentHeat_evap);
    } else if (mode == HEATING_MODE) {
        Qc = calculateHeatAbsorbed(h1_adjusted, h4_adjusted, m, latentHeat_evap);
        Qe = calculateHeatRejected(h2_adjusted, h3_adjusted, m, latentHeat_cond);
    }

    // Work done by the expansion valve
    float We = calculateWorkDone(h3_adjusted, h4_adjusted, m);

    // Coefficient of Performance (COP)
    float COP = calculateCOP(Qe, Wc_actual);

    // Second-law efficiencies
    float secondLawEfficiency_compressor = (Wc_ideal > 0) ? calculateSecondLawEfficiency(Wc_actual, Wc_ideal) : 0.0;
    float secondLawEfficiency_turbine = (Wt_ideal > 0) ? calculateSecondLawEfficiency(Wt_actual, Wt_ideal) : 0.0;

    // Heat transfer coefficients with fouling factors
    float Re_evap = (rho * v * D) / 1.0; // Reynolds number in evaporator
    float Re_cond = (rho * v * D) / 1.0; // Reynolds number in condenser
    float Pr_evap = 7.0; // Prandtl number in evaporator
    float Pr_cond = 7.0; // Prandtl number in condenser

    float h_evap = calculateHeatTransferCoefficient(Re_evap, Pr_evap, D, foulingFactor_evap);
    float h_cond = calculateHeatTransferCoefficient(Re_cond, Pr_cond, D, foulingFactor_cond);

    float A_evap = 5.0; // Area in evaporator (m²); assume a value
    float A_cond = 5.0; // Area in condenser (m²); assume a value
    float dT_evap = T_evap - T_cond; // Temperature difference in evaporator (°C)
    float dT_cond = T_cond - T_evap; // Temperature difference in condenser (°C)

    float Q_evap = calculateHeatTransferRate(h_evap, A_evap, dT_evap);
    float Q_cond = calculateHeatTransferRate(h_cond, A_cond, dT_cond);

    float P_v = 1.2; // Partial pressure of water vapor (bar) 
    float RH = 60.0; // Relative humidity (%)
    float P = 1.0; // Total atmospheric pressure (bar)

    float humidityRatio = calculateHumidityRatio(P_v, P);
    float dewPoint = calculateDewPoint(T_evap, RH);
    float wetBulbTemperature = calculateWetBulbTemperature(T_evap, dewPoint);

    printf("Mechanical Refrigeration Cycle");
    printf("\n------------------------------\n");
    printf("\nHeat Pump Mode: %s\n", (mode == COOLING_MODE) ? "Cooling" : "Heating");
    printf("\nState 1 (Evaporator Inlet - Saturated Liquid):\n");
    printf("  Temperature = %.2f°C\n", T_evap);
    printf("  Pressure    = %.2f bar\n", P_evap);
    printf("  Enthalpy    = %.2f kJ/kg\n", h1_liquid);
    printf("  Cp = %.2f kJ/(kg·K)\n", getSpecificHeatCapacity(T_evap, P_evap, 1));
    printf("  Cv = %.2f kJ/(kg·K)\n", getSpecificHeatCapacity(T_evap, P_evap, 0));
    printf("\nState 2 (Compressor Outlet - Saturated Vapor):\n");
    printf("  Temperature = %.2f°C\n", T_cond);
    printf("  Pressure    = %.2f bar\n", P_cond);
    printf("  Enthalpy    = %.2f kJ/kg\n", h2_vapor);
    printf("  Cp = %.2f kJ/(kg·K)\n", getSpecificHeatCapacity(T_cond, P_cond, 1));
    printf("  Cv = %.2f kJ/(kg·K)\n", getSpecificHeatCapacity(T_cond, P_cond, 0));
    printf("\nState 3 (Condenser Outlet - Saturated Liquid):\n");
    printf("  Temperature = %.2f°C\n", T_cond);
    printf("  Pressure    = %.2f bar\n", P_cond);
    printf("  Enthalpy    = %.2f kJ/kg\n", h3_liquid);
    printf("  Cp = %.2f kJ/(kg·K)\n", getSpecificHeatCapacity(T_cond, P_cond, 1));
    printf("  Cv = %.2f kJ/(kg·K)\n", getSpecificHeatCapacity(T_cond, P_cond, 0));
    printf("\nState 4 (Expansion Valve Outlet - Saturated Vapor):\n");
    printf("  Temperature = %.2f°C\n", T_evap);
    printf("  Pressure    = %.2f bar\n", P_evap);
    printf("  Enthalpy    = %.2f kJ/kg\n", h4_vapor);
    printf("  Cp = %.2f kJ/(kg·K)\n", getSpecificHeatCapacity(T_evap, P_evap, 1));
    printf("  Cv = %.2f kJ/(kg·K)\n", getSpecificHeatCapacity(T_evap, P_evap, 0));
    printf("\nPressure Drops:\n");
    printf("  Evaporator Pressure Drop   = %.2f kPa\n", deltaP_e);
    printf("  Condenser Pressure Drop    = %.2f kPa\n", deltaP_c);
    printf("  Expansion Valve Pressure Drop = %.2f kPa\n", deltaP_v);
    printf("\nAdjusted Enthalpies:\n");
    printf("  State 1 (Adjusted)         = %.2f kJ/kg\n", h1_adjusted);
    printf("  State 2 (Adjusted)         = %.2f kJ/kg\n", h2_adjusted);
    printf("  State 3 (Adjusted)         = %.2f kJ/kg\n", h3_adjusted);
    printf("  State 4 (Adjusted)         = %.2f kJ/kg\n", h4_adjusted);
    printf("\nHeat Transfer Coefficients:\n");
    printf("  Evaporator Heat Transfer Coefficient (h_evap) = %.2f W/(m²·K)\n", h_evap);
    printf("  Condenser Heat Transfer Coefficient (h_cond) = %.2f W/(m²·K)\n", h_cond);
    printf("\nHeat Transfer Rates:\n");
    printf("  Heat Transfer Rate in Evaporator (Q_evap) = %.2f kJ/s\n", Q_evap / 1000); // Convert W to kJ/s
    printf("  Heat Transfer Rate in Condenser (Q_cond) = %.2f kJ/s\n", Q_cond / 1000); // Convert W to kJ/s
    printf("\nPerformance Metrics:\n");
    printf("  Ideal work done by compressor (Wc_ideal): %.2f kJ/s\n", Wc_ideal);
    printf("  Actual work done by compressor (Wc_actual): %.2f kJ/s\n", Wc_actual);
    printf("  Second-law efficiency of compressor: %.2f\n", secondLawEfficiency_compressor);
    printf("  Ideal work done by turbine (Wt_ideal): %.2f kJ/s\n", Wt_ideal);
    printf("  Actual work done by turbine (Wt_actual): %.2f kJ/s\n", Wt_actual);
    printf("  Second-law efficiency of turbine: %.2f\n", secondLawEfficiency_turbine);
    printf("  Heat rejected in condenser (Qc): %.2f kJ/s\n", Qc);
    printf("  Heat absorbed in evaporator (Qe): %.2f kJ/s\n", Qe);
    printf("  Work done by expansion valve (We): %.2f kJ/s\n", We);
    printf("  Coefficient of Performance (COP): %.2f\n", COP);
    printf("\nHumidity Calculations:\n");
    printf("  Humidity Ratio (ω): %.2f kg/kg\n", humidityRatio);
    printf("  Dew Point Temperature (T_dp): %.2f°C\n", dewPoint);
    printf("  Wet Bulb Temperature (T_wb): %.2f°C\n", wetBulbTemperature);

    return 0;
}

float calculateWorkDone(float h_out, float h_in, float m) {
    return m * (h_out - h_in);
}

float calculateHeatRejected(float h_out, float h_in, float m, float latentHeat) {
    return m * (h_out - h_in) + latentHeat; 
}

float calculateHeatAbsorbed(float h_out, float h_in, float m, float latentHeat) {
    return m * (h_out - h_in) - latentHeat; 
}

float calculateCOP(float Qe, float Wc_actual) {
    return Qe / Wc_actual;
}

float calculateActualCompressorWork(float h2s, float h1, float efficiency, float m) {
    float Wc_ideal = calculateWorkDone(h2s, h1, m);
    return (Wc_ideal > 0) ? Wc_ideal / efficiency : 0.0; 
}

float calculateActualTurbineWork(float h2, float h1s, float efficiency, float m) {
    float Wt_ideal = calculateWorkDone(h1s, h2, m);
    return (Wt_ideal > 0) ? Wt_ideal * efficiency : 0.0; 
}

float adjustEnthalpyForPressureDrop(float h, float deltaP) {
    return h - (deltaP * 0.01); 
}

float getSaturatedLiquidEnthalpy(float T, float P) {
    if (T == 5.0 && P == 1.2) return 250.0; 
    if (T == 35.0 && P == 4.5) return 300.0; 
    return 0.0;
}

float getSaturatedVaporEnthalpy(float T, float P) {
    if (T == 5.0 && P == 1.2) return 450.0; 
    if (T == 35.0 && P == 4.5) return 500.0; 
    return 0.0;
}

float getLatentHeat(float T, float P) {
    if (T == 5.0 && P == 1.2) return 200.0; 
    if (T == 35.0 && P == 4.5) return 220.0; 
    return 0.0;
}

float calculatePressureDrop(float L, float D, float rho, float v, float f) {
    // Darcy-Weisbach equation for pressure drop
    return f * (L / D) * (0.5 * rho * v * v) / 1000; // Pressure drop in kPa
}

float calculateExpansionPressureDrop(float P1, float T1, float T2) {
    // Model for expansion valve pressure drop
    return (P1 - (T2 / T1) * P1) * 100; // Convert to kPa
}

float getSpecificHeatCapacity(float T, float P, int isConstantPressure) {
    if (isConstantPressure) return 4.18; 
    else return 3.98; 
}

float calculateIsentropicEnthalpy(float h1, float P1, float P2) {
    return h1 * (P2 / P1); 
}

float calculateIdealWork(float h1, float h2, float m) {
    return m * (h2 - h1);
}

float calculateSecondLawEfficiency(float W_actual, float W_ideal) {
    if (W_ideal == 0) return 0.0;
    return W_actual / W_ideal;
}

float calculateHeatTransferCoefficient(float Re, float Pr, float D, float foulingFactor) {
    // Dittus-Boelter equation for turbulent flow
    float h = 0.023 * pow(Re, 0.8) * pow(Pr, 0.3) * (THERMAL_CONDUCTIVITY_WATER / D); 
    return h * (1 + foulingFactor); 
}

float calculateHeatTransferRate(float h, float area, float dT) {
    return h * area * dT;
}

float calculateHumidityRatio(float P_v, float P) {
    return 0.622 * (P_v / (P - P_v));
}

float calculateDewPoint(float T, float RH) {
    float alpha = (DEW_POINT_A * T) / (DEW_POINT_B + T) + log(RH / 100.0);
    return (DEW_POINT_B * alpha) / (DEW_POINT_A - alpha);
}

float calculateWetBulbTemperature(float T, float T_dp) {
    return T - (T - T_dp) / 3.3;
}
