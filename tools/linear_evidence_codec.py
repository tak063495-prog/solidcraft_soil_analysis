"""Lossless float64 storage of actual dense solve records, no quantization."""
import gzip
import json
import struct
import numpy as np

MAGIC = b'SOIL-LINEAR-F64-1\n'
ARRAYS = ('matrix', 'rhs', 'solution')


def encode_rows(rows, destination):
    previous = {}
    with destination.open('xb') as raw, gzip.GzipFile(fileobj=raw, filename='', mode='wb', compresslevel=6, mtime=0) as out:
        out.write(MAGIC)
        for row in rows:
            arrays = {key: np.asarray(row[key], dtype='<f8') for key in ARRAYS}
            if any(not np.isfinite(a).all() for a in arrays.values()):
                raise ValueError('nonfinite actual linear system')
            header = dict(keys=list(row), metadata={k: v for k, v in row.items() if k not in ARRAYS},
                          shapes={k: a.shape for k, a in arrays.items()}, xor_previous=True)
            b = json.dumps(header, allow_nan=False).encode()
            out.write(struct.pack('<I', len(b))); out.write(b)
            for key, a in arrays.items():
                bits = a.view('<u8')
                old = previous.get(key)
                payload = np.bitwise_xor(bits, old) if old is not None and old.shape == bits.shape else bits
                out.write(payload.tobytes(order='C'))
                previous[key] = bits.copy()


def encode(source, destination):
    with gzip.open(source, 'rt', encoding='utf-8') as src:
        encode_rows((json.loads(line) for line in src), destination)


def read_rows(path):
    def exact(stream, size):
        data = stream.read(size)
        if len(data) != size:
            raise ValueError('truncated lossless linear evidence')
        return data
    with gzip.open(path, 'rb') as stream:
        previous = {}
        if exact(stream, len(MAGIC)) != MAGIC:
            raise ValueError('unknown linear evidence encoding')
        while True:
            prefix = stream.read(4)
            if not prefix:
                return
            if len(prefix) != 4:
                raise ValueError('truncated linear evidence header')
            size, = struct.unpack('<I', prefix)
            if size > 1024*1024:
                raise ValueError('oversized linear evidence header')
            header = json.loads(exact(stream, size)); arrays = {}
            for key in ARRAYS:
                shape = header['shapes'][key]
                if not shape or len(shape) > 2 or any(type(n) is not int or not 0 < n <= 10000 for n in shape):
                    raise ValueError('invalid actual linear evidence shape')
                count = int(np.prod(shape))
                bits = np.frombuffer(exact(stream, count*8), dtype='<u8').reshape(shape)
                old = previous.get(key)
                if header.get('xor_previous') and old is not None and old.shape == bits.shape:
                    bits = np.bitwise_xor(bits, old)
                previous[key] = bits.copy()
                arrays[key] = bits.view('<f8').tolist()
            yield {key: arrays[key] if key in ARRAYS else header['metadata'][key] for key in header['keys']}


def json_lines(path):
    for row in read_rows(path):
        yield (json.dumps(row, allow_nan=False)+'\n').encode()
