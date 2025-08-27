#!/bin/bash
# encode_analysis/scripts/setup_workspace.sh
#
# Purpose: Set up a robust analysis environment for variant calling pipelines.
# This script installs all necessary tools, configures persistent workspace settings,
# and handles network interruptions and long-running processes.
#
# Note: Run this script as root.

# Exit immediately if a command exits with a non-zero status, and print commands as they execute
set -euo pipefail
set -x

# Check for root privileges
if [[ "$EUID" -ne 0 ]]; then
  echo "Please run this script as root."
  exit 1
fi

echo "Setting up analysis environment..."

# Base directory for the workspace
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis"

# Create comprehensive directory structure
mkdir -p "${BASE_DIR}"/{scripts,logs,configs,pipeline/tools}
mkdir -p "${BASE_DIR}/data/variants"/{deepsomatic,deepvariant,cell_line_analysis_v{3..5}}
mkdir -p "${BASE_DIR}/data"/{bam,fastq,reference}

# Update system package list and install essential tools
echo "Updating system and installing essential tools..."
apt-get update
apt-get install -y \
    build-essential \
    wget \
    curl \
    git \
    screen \
    tmux \
    htop \
    vim \
    less \
    parallel \
    psmisc \
    openjdk-11-jdk \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev

# Install Python packages via apt-get
echo "Installing system-level Python packages..."
apt-get install -y \
    python3-pip \
    python3-yaml \
    python3-psutil \
    python3-tqdm \
    python3-numpy \
    python3-pandas \
    python3-matplotlib

# Upgrade pip (optional but recommended) and install additional Python packages via pip
echo "Installing Python packages via pip..."
python3 -m pip install --upgrade pip
python3 -m pip install --no-cache-dir \
    pyyaml \
    psutil \
    tqdm \
    requests \
    pandas \
    matplotlib \
    seaborn \
    biopython

# Install bioinformatics tools
echo "Installing bioinformatics tools..."
apt-get install -y \
    bcftools \
    samtools \
    tabix \
    bwa \
    minimap2 \
    bedtools \
    vcftools

# Set up Picard tools
echo "Setting up Picard tools..."
PICARD_VERSION="2.27.5"
PICARD_JAR="${BASE_DIR}/pipeline/tools/picard.jar"

if [ ! -f "${PICARD_JAR}" ]; then
    wget -O "${PICARD_JAR}" "https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar"
fi

# Create Picard wrapper script
cat > /usr/local/bin/picard << 'EOF'
#!/bin/bash
java -Xmx8g -jar /mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/pipeline/tools/picard.jar "$@"
EOF
chmod +x /usr/local/bin/picard

# Create symlink for the yaml package if needed (only if the file doesn't exist)
if [ ! -e /usr/lib/python3/dist-packages/yaml ]; then
    ln -s /usr/local/lib/python3/dist-packages/yaml /usr/lib/python3/dist-packages/yaml
fi

# Set up environment variables
echo "Setting up environment variables..."
cat > /etc/profile.d/encode_analysis.sh << 'EOF'
export ENCODE_ANALYSIS_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis"
export PATH="$ENCODE_ANALYSIS_DIR/scripts:$PATH"
export PYTHONPATH="$ENCODE_ANALYSIS_DIR:$PYTHONPATH"
export JAVA_OPTS="-Xmx16g"
EOF
chmod +x /etc/profile.d/encode_analysis.sh

# Set up screen configuration for persistent sessions
echo "Configuring screen for persistent sessions..."
cat > /root/.screenrc << 'EOF'
# Enable scroll buffer
defscrollback 10000

# Status line at bottom
hardstatus alwayslastline
hardstatus string '%{= kG}[ %{G}%H %{g}][%= %{= kw}%?%-Lw%?%{r}(%{W}%n*%f%t%?(%u)%?%{r})%{w}%?%+Lw%?%?%= %{g}][%{B} %m-%d %{W}%c %{g}]'

# Enable mouse scrolling
termcapinfo xterm* ti@:te@

# Don't display the copyright page
startup_message off

# Allow bold colors
attrcolor b ".I"

# Tell screen how to set colors
termcapinfo xterm 'Co#256:AB=\E[48;5;%dm:AF=\E[38;5;%dm'

# Enable 256 color terminal
term screen-256color

# Cache 10000 lines for scroll back
defscrollback 10000

# Set up logging
logfile /var/log/screen/screenlog-%t-%n-%Y%m%d-%c
logtstamp on
logtstamp after 5
EOF

# Create screen log directory
mkdir -p /var/log/screen
chmod 777 /var/log/screen

# Create directories for individual cell line analyses
for cell_line in K562 HepG2 Caki2 A549 GM23248 NCI-H460; do
    echo "Setting up directory structure for ${cell_line}..."
    mkdir -p "${BASE_DIR}/data/variants/cell_line_analysis_v5/${cell_line}"/{intermediate,final}
    mkdir -p "${BASE_DIR}/data/bam/${cell_line}"
done
