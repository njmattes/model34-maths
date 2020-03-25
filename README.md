# Obstructures Model 34 maths.

This repo contains various maths required for the design of Obstructures'
Model 34 bass guitar.

If you want to run this and have no idea what you're doing, open
Terminal (assuming you're on OS X) and do the following:

1. Clone this repo onto your machine
   ```
   git clone https://github.com/njmattes/model34-maths.git
   cd model34-maths
   ```

2. Install the required python packages
   ```
   pip install -r requirement.txt
   ```

3. Add the directory to your python path and execute the main.py file
   ```
   PYTHONPATH="$PYTHONPATH:./" python model34_maths/main.py
   ```

As the optimizer runs, you'll see values for the objective function appear.
If the optimizer is able to successfully find an optimal solution, you'll
eventually see areas and deflection for the original neck and the optimal
neck appear, followed by optimal parameters.

