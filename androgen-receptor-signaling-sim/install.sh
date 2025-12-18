#!/bin/bash

# Androgen Receptor Signaling Simulation - Installation and Test Script

echo "=========================================="
echo "AR SIGNALING SIMULATION - SETUP"
echo "=========================================="
echo ""

# Check Python version
echo "1. Checking Python version..."
python3 --version

# Create virtual environment (optional)
echo ""
echo "2. Creating virtual environment (optional)..."
python3 -m venv venv 2>/dev/null || echo "Virtual environment already exists or creation skipped"

# Activate virtual environment
echo ""
echo "3. To activate virtual environment, run:"
echo "   source venv/bin/activate  # on macOS/Linux"
echo "   venv\\Scripts\\activate    # on Windows"
echo ""
read -p "Press enter to continue..."

# Install dependencies
echo ""
echo "4. Installing dependencies..."
pip install -r requirements.txt

# Install package
echo ""
echo "5. Installing package..."
pip install -e .

# Create output directories
echo ""
echo "6. Creating output directories..."
mkdir -p plots results data/processed

# Run basic test
echo ""
echo "7. Running basic test simulation..."
python scripts/run_simulation.py

echo ""
echo "=========================================="
echo "SETUP COMPLETE!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "  1. Check plots/ directory for generated figures"
echo "  2. Run: python scripts/drug_comparison.py"
echo "  3. Run: python scripts/sensitivity_analysis.py"
echo "  4. Launch dashboard: python -m src.visualization.interactive"
echo "  5. Explore notebooks: jupyter notebook notebooks/"
echo ""
echo "Documentation:"
echo "  - readme.md - comprehensive guide"
echo "  - QUICKSTART.md - quick reference"
echo "  - PROJECT_SUMMARY.md - technical overview"
echo ""
