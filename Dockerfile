FROM python:3.8

SHELL ["/bin/sh", "-c"]

# install OS packages
RUN apt-get update && \
    apt-get install -yq git wget gdebi-core build-essential software-properties-common \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt-get update -y && apt-get install google-cloud-cli -y

# install CellphoneDB
RUN pip3 install git+https://github.com/dpeerlab/spectra.git --no-cache-dir
RUN pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

RUN ln -s /usr/bin/python3 /usr/bin/python