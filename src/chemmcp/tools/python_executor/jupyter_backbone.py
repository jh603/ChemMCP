import os
import time
import uuid
import json
import logging
import tempfile
import requests
import atexit
import signal
import re
from websocket import create_connection, WebSocketTimeoutException

logger = logging.getLogger(__name__)


DOCKERFILE = """\
FROM continuumio/miniconda3:latest
RUN apt-get update && apt-get install -y \
    cmake \
    pkg-config \
    build-essential \
    && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir jupyter jupyter_kernel_gateway lmdb>=1.6.2 molbloom>=2.3.4 pandas>=2.2.3 pubchempy>=1.0.4 rdchiral>=1.1.0 rdkit>=2024.9.6 requests>=2.32.3 rxn4chemistry>=1.14.0 scikit-learn>=1.6.1 selfies>=2.2.0 sentencepiece>=0.2.0 tavily-python>=0.7.0 transformers>=4.51.3 matplotlib
"""


def strip_ansi(text: str) -> str:
    """Remove ANSI escape sequences for clean output."""
    ansi_re = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')
    return ansi_re.sub('', text)

class KernelSession:
    """
    Minimal Jupyter kernel session over HTTP/WebSocket, no external dependencies.
    """
    def __init__(self, base_url: str, convid: str, kernel_name: str = "python3"):
        self.base_url = base_url.rstrip('/')
        self.ws_url = self.base_url.replace('http://', 'ws://').replace('https://', 'wss://')
        self.convid = convid

        # Connectivity check
        try:
            resp = requests.get(f"{self.base_url}/api", timeout=5)
            resp.raise_for_status()
        except requests.RequestException as e:
            raise ConnectionError(
                f"Cannot reach Jupyter Kernel Gateway at {self.base_url}. "
                "Ensure the gateway is running and the port is exposed (e.g., `docker run -p 8888:8888 ...`). "
                f"Error: {e}")

        # Start kernel
        self.kernel_id = self._start_kernel(kernel_name)
        # Connect websocket
        self.ws = self._connect_ws()
        # Initial setup
        self._setup()

    def _start_kernel(self, kernel_name: str) -> str:
        url = f"{self.base_url}/api/kernels"
        resp = requests.post(url, json={"name": kernel_name})
        resp.raise_for_status()
        kid = resp.json()["id"]
        logger.info(f"Started kernel {kid}")
        return kid

    def _connect_ws(self):
        ws_endpoint = f"{self.ws_url}/api/kernels/{self.kernel_id}/channels"
        for _ in range(10):
            try:
                ws = create_connection(ws_endpoint, timeout=5)
                logger.info("WebSocket connected to kernel channels")
                return ws
            except Exception:
                time.sleep(0.5)
        raise RuntimeError("Failed to open WebSocket to kernel channels")

    def _setup(self):
        # Disable ANSI coloring for consistency
        self.execute("%colors nocolor")

    def execute(self, code: str, timeout: int = 60) -> str:
        msg_id = uuid.uuid4().hex
        header = {"msg_id": msg_id, "username": "", "session": self.convid, "msg_type": "execute_request", "version": "5.0"}
        message = {
            "header": header,
            "parent_header": {},
            "metadata": {},
            "content": {"code": code, "silent": False, "store_history": False, "user_expressions": {}, "allow_stdin": False},
            "channel": "shell"
        }
        self.ws.send(json.dumps(message))

        output = []
        start = time.time()
        while True:
            if time.time() - start > timeout:
                raise TimeoutError(f"Execution timed out after {timeout}s")
            try:
                raw = self.ws.recv()
            except WebSocketTimeoutException:
                continue
            msg = json.loads(raw)
            if msg.get('parent_header', {}).get('msg_id') != msg_id:
                continue
            mtype = msg.get('msg_type')
            content = msg.get('content', {})
            if mtype == 'stream':
                output.append(content.get('text', ''))
            elif mtype in ('execute_result', 'display_data'):
                data = content.get('data', {}).get('text/plain')
                if data:
                    output.append(data)
            elif mtype == 'error':
                tb = '\n'.join(content.get('traceback', []))
                output.append(tb)
                break
            elif mtype == 'execute_reply':
                break

        result = strip_ansi(''.join(output)).strip()
        return result or '[No output]'

    def shutdown(self):
        # Stop kernel and close WS
        try:
            requests.delete(f"{self.base_url}/api/kernels/{self.kernel_id}")
            logger.info(f"Shutdown kernel {self.kernel_id}")
        except Exception as e:
            logger.debug(f"Error shutting down kernel: {e}")
        try:
            self.ws.close()
        except Exception:
            pass

class JupyterBackbone:
    """
    Runs code in a Docker-hosted or existing Jupyter kernel gateway, auto-cleaning containers on exit.

    :param kernel_url: connect to this gateway if given
    :param image: Docker image name
    :param mem_limit: memory limit
    :param cpus: CPU cores limit
    """
    DEFAULT_IMAGE = "chemmcp-python-executor:latest"
    JUPYTER_PORT = 8888

    def __init__(self, kernel_url: str = None, image: str = None, mem_limit: str = "8g", cpus: float = 2.0):
        self.conv_id = str(uuid.uuid4())
        self.image = image or self.DEFAULT_IMAGE
        self.mem_limit = mem_limit
        self.cpus = cpus
        self.container = None
        self.session = None

        # Register cleanup on normal exit and signals
        atexit.register(self.close)
        for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
            try:
                signal.signal(sig, lambda s, f: self.close())
            except Exception:
                pass

        if kernel_url:
            self.kernel_url = kernel_url.rstrip('/')
            logger.info(f"Connecting to gateway at {self.kernel_url}")
            self.session = KernelSession(self.kernel_url, self.conv_id)
        else:
            import docker
            client = docker.from_env()
            # Ensure or build image
            try:
                client.images.get(self.image)
                logger.info(f"Docker image {self.image} found")
            except docker.errors.ImageNotFound:
                logger.info(f"Building Docker image {self.image}")
                tmp = tempfile.mkdtemp()
                with open(os.path.join(tmp, 'Dockerfile'), 'w') as df:
                    df.write(DOCKERFILE)
                client.images.build(path=tmp, tag=self.image)
                logger.info("Image built locally")
            # Launch container
            ports = {f"{self.JUPYTER_PORT}/tcp": None}
            self.container = client.containers.run(
                self.image,
                command=f"jupyter kernelgateway --ip=0.0.0.0 --port={self.JUPYTER_PORT}",
                detach=True,
                remove=True,
                ports=ports,
                mem_limit=self.mem_limit,
                nano_cpus=int(self.cpus * 1e9)
            )
            self.container.reload()
            binding = self.container.attrs['NetworkSettings']['Ports'][f"{self.JUPYTER_PORT}/tcp"][0]
            port = int(binding['HostPort'])
            self.kernel_url = f"http://localhost:{port}"
            logger.info(f"Gateway container started at {self.kernel_url}")
            time.sleep(2)
            self.session = KernelSession(self.kernel_url, self.conv_id)

    def run_code(self, code: str, timeout: int = 600) -> str:
        logger.info(f"Executing code (conv_id={self.conv_id}): {code}")
        out = self.session.execute(code, timeout)
        logger.info(f"Result: {out}")
        return out

    def close(self):
        if self.session:
            self.session.shutdown()
            self.session = None
        if self.container:
            try:
                logger.info(f"Stopping container {self.container.id}")
                self.container.stop()
            except Exception as e:
                logger.debug(f"Error stopping container: {e}")
            self.container = None

    def __del__(self):
        self.close()
