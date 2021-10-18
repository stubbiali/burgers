#!/bin/bash

PYTHON=python
VENV=venv
PIP_UPGRADE=1
FRESH_INSTALL=0

function install()
{
  # activate environment
  source $VENV/bin/activate

  # upgrade pip and setuptools
  if [ "$PIP_UPGRADE" -gt 0 ]; then
    pip install --upgrade pip setuptools wheel
  fi

  # install requirements
  cat requirements.txt | xargs -n 1 pip install

  # install gt sources
  python -m gt4py.gt_src_manager install

  # deactivate environment
  deactivate

  # On OSX only: change matplotlib backend from macosx to TkAgg
  if [[ "$OSTYPE" == "darwin"* ]]; then
    cat $VENV/lib/$PYTHON/site-packages/matplotlib/mpl-data/matplotlibrc | \
      sed -e 's/^backend.*: macosx/backend : TkAgg/g' > /tmp/.matplotlibrc && \
      cp /tmp/.matplotlibrc $VENV/lib/$PYTHON/site-packages/matplotlib/mpl-data/matplotlibrc && \
      rm /tmp/.matplotlibrc
  fi
}

if [ "$FRESH_INSTALL" -gt 0 ]
then
  echo -e "Creating new environment..."
  rm -rf $VENV
  $PYTHON -m venv $VENV
fi

install || deactivate

echo -e ""
echo -e "Command to activate environment:"
echo -e "\t\$ source $VENV/bin/activate"
echo -e ""
echo -e "Command to deactivate environment:"
echo -e "\t\$ deactivate"
echo -e ""
