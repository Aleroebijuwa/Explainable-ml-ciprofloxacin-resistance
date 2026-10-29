#!/usr/bin/env bash
# Installs everything needed by extract_features.py on a fresh Ubuntu 22.04/24.04 machine (e.g. AWS EC2).
# Usage: bash aws/setup.sh
set -euo pipefail

echo "Installing system packages..."
sudo apt-get update -qq
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
    bwa samtools curl git make gcc zlib1g-dev python3-venv python3-pip > /dev/null

echo "Installing Python packages..."
python3 -m venv "$HOME/venv"
"$HOME/venv/bin/pip" install -q --upgrade pip
"$HOME/venv/bin/pip" install -q pandas numpy pysam biopython requests

if ! command -v kma > /dev/null; then
    echo "Building KMA..."
    rm -rf /tmp/kma
    git clone -q --depth 1 https://bitbucket.org/genomicepidemiology/kma.git /tmp/kma
    make -s -C /tmp/kma > /dev/null 2>&1
    sudo cp /tmp/kma/kma /tmp/kma/kma_index /usr/local/bin/
fi

if ! command -v mash > /dev/null; then
    echo "Installing Mash..."
    curl -sL -o /tmp/mash.tar https://github.com/marbl/Mash/releases/download/v2.3/mash-Linux64-v2.3.tar
    tar -xf /tmp/mash.tar -C /tmp
    sudo cp /tmp/mash-Linux64-v2.3/mash /usr/local/bin/
fi

echo
echo "Tool versions (record these in the methods section):"
bwa 2>&1 | grep Version || true
samtools --version | head -1
echo "kma $(kma -v 2>&1 | head -1)"
echo "mash $(mash --version)"
echo "CPUs: $(nproc)   Free disk: $(df -h --output=avail "$HOME" | tail -1)"
echo
echo "Setup complete."
