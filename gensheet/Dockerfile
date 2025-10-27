FROM archlinux:latest

# Install dependencies
RUN pacman -Syu --noconfirm

# Install Pixi
RUN pacman -S --noconfirm pixi

# Copy pixi.lock and pixi.toml
COPY pixi.lock pixi.toml ./

# Run Pixi to install dependencies
RUN pixi install

# Expose the port Pixi uses
EXPOSE 8501

# Activate the Pixi environment
RUN pixi shell-hook > /shell-hook.sh
RUN echo 'exec "$@"' >> /shell-hook.sh
COPY app.py app.py
# Set the entrypoint to run Pixi
ENTRYPOINT ["/bin/bash", "/shell-hook.sh", "streamlit", "run", "app.py", "--server.maxUploadSize", "200000"]
#ENTRYPOINT [ "pixi run streamlit run app.py --server.maxUploadSize 200000" ]