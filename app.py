from flask import Flask, render_template, request, jsonify
import numpy as np
import matplotlib.pyplot as plt
from set_stack import set_stack
from ATR1D import ATR1D

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/calculate', methods=['POST'])
def calculate():
    data = request.json
    
    # Extract data from the request
    lam = np.array(data.get('lam', np.linspace(380, 1080, 851)))
    Bragg2 = data.get('bragg1', ["SIO2.csv","ZrO2.csv","SIO2.csv","ZrO2.csv","SIO2.csv","ZrO2.csv","SIO2.csv","ZrO2.csv"])
    dBragg2 = data.get('dBragg1', [106.36,29.51,13.06,99.32,16.04,30.71,44.67,11.84])
    Bragg1 = data.get('bragg2')
    dBragg1 = data.get('dBragg2')
    dgls = data.get('dgls', 320000)
    theta = data.get('theta', 0)

    # General input data
    matPSC = Bragg1 + ["GLS_NEW.csv"]
    dPSC = dBragg1 + [dgls]

    # Add Bragg2 if it exists
    if Bragg2 is not None and dBragg2 is not None:
        matPSC += Bragg2
        dPSC += dBragg2

    incoh = 1e3  # incoherent layer is 1000nm thickness
    dAZO = 0
    dEVA = 0
    dARC = 0

    # PSC with AZO and EVA on PV
    materials = ["air"] + matPSC + ["air", "air", "air", "air"]
    d = dPSC + [dAZO, dEVA, dARC]

    try:
        stack = set_stack(materials, d, lam, theta, incoh)
        _, T_PSC_AZO_EVA_PV, _ = ATR1D(stack)
        print("Calculation completed successfully")
    except Exception as e:
        print("Error during calculation:", str(e))
        return jsonify({"error": str(e)}), 500

    # Prepare response
    response = {
        'wavelengths': lam.tolist(),
        'transmittance': T_PSC_AZO_EVA_PV['sp'].tolist()
    }

    return jsonify(response)

if __name__ == '__main__':
    app.run(debug=True)
