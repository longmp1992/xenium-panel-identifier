# Xenium Panel Identifier - image for SciLifeLab Serve ("Gradio app", port 7860, user id 1000)
FROM python:3.12-slim

ENV USER=appuser \
    HOME=/home/appuser \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1
RUN useradd -m -u 1000 $USER

WORKDIR $HOME/app

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY main.py panel_core.py extract_features.py start-script.sh ./
COPY data ./data
RUN chmod +x start-script.sh && mkdir -p $HOME/app/tmp && chown -R $USER:$USER $HOME

USER $USER

# Serve limits uploads to 100 MB; the app shows a clear message above that size.
# Override with -e MAX_UPLOAD_MB=... (or an empty value for no app-side limit) when running elsewhere.
ENV GRADIO_SERVER_NAME=0.0.0.0 \
    GRADIO_SERVER_PORT=7860 \
    GRADIO_TEMP_DIR=/home/appuser/app/tmp \
    GRADIO_ANALYTICS_ENABLED=False \
    MAX_UPLOAD_MB=100

EXPOSE 7860
CMD ["./start-script.sh"]
