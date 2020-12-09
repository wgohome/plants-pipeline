# Downloads kallisto if not present in ./dependencies
if [ ! -d dependencies/kallisto ]; then
  echo "Downloading kallisto ..."
  # Specific download url for mac:
  # wget -O - https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_mac-v0.46.1.tar.gz | tar -C dependencies/ -xf -
  # For linux:
  wget -O kallisto.tar.gz https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz
  tar -C dependencies/ -xf kallisto.tar.gz
  rm -r kallisto.tar.gz
else
  echo "kallisto already present in ./dependencies"
fi

# Downloads aspera CLI if not present in native applications dir
if [ ! -d ~/.aspera ]; then
  echo "Downloading Aspera CLI ..."
  # For mac
  # wget -O aspera_cli.sh https://download.asperasoft.com/download/sw/cli/3.9.2/ibm-aspera-cli-3.9.2.1426.c59787a-mac-10.7-64-release.sh
  # For linux
  wget -O aspera_cli.sh https://download.asperasoft.com/download/sw/cli/3.9.2/ibm-aspera-cli-3.9.2.1426.c59787a-linux-64-release.sh
  sh aspera_cli.sh
  rm aspera_cli.sh
else
  echo "aspera CLI already present in ~/.aspera"
fi

PATH=`pwd`/dependencies/kallisto:$PATH
PATH=~/.aspera/cli/bin:$PATH # Specific path to linux
export $PATH
echo "\$PATH exported, 'kallisto' and 'ascp' can be called from the CLI"

# Install chromedriver
# sudo apt-get install chromium-chromedriver
