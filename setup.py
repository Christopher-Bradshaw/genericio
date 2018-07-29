from distutils.core import setup, Extension
from Cython.Build import cythonize

setup(
        # exts to build - hence "build_ext" in the compile arg
        ext_modules = cythonize([
            Extension(
                "wrapper",
                sources=["wrapper.pyx"],
                language="c++",
                extra_link_args=["-fopenmp"],
                extra_compile_args=["-DGENERICIO_NO_MPI"],
                extra_objects=[
                    "./frontend/GenericIO.o",
                    "./frontend/thirdparty/blosc/blosc.o",
                    "./frontend/thirdparty/blosc/blosclz.o",
                    "./frontend/thirdparty/blosc/shuffle.o",
                    "./frontend/thirdparty/blosc/bitshuffle-generic.o",
                    "./frontend/thirdparty/blosc/shuffle-generic.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zlib-1.2.8/gzwrite.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zlib-1.2.8/crc32.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zlib-1.2.8/inffast.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zlib-1.2.8/zutil.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zlib-1.2.8/infback.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zlib-1.2.8/deflate.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zlib-1.2.8/inflate.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zlib-1.2.8/gzread.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zlib-1.2.8/gzlib.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zlib-1.2.8/gzclose.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zlib-1.2.8/uncompr.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zlib-1.2.8/compress.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zlib-1.2.8/inftrees.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zlib-1.2.8/trees.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zlib-1.2.8/adler32.o",
                    "./frontend/thirdparty/blosc/internal-complibs/lz4-1.7.2/lz4.o",
                    "./frontend/thirdparty/blosc/internal-complibs/lz4-1.7.2/lz4hc.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zstd-0.7.4/legacy/zstd_v01.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zstd-0.7.4/legacy/zstd_v02.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zstd-0.7.4/legacy/zstd_v03.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zstd-0.7.4/legacy/zstd_v06.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zstd-0.7.4/legacy/zstd_v04.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zstd-0.7.4/legacy/zstd_v05.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zstd-0.7.4/compress/fse_compress.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zstd-0.7.4/compress/zstd_compress.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zstd-0.7.4/compress/huf_compress.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zstd-0.7.4/compress/zbuff_compress.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zstd-0.7.4/common/entropy_common.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zstd-0.7.4/common/xxhash.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zstd-0.7.4/common/zstd_common.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zstd-0.7.4/common/fse_decompress.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zstd-0.7.4/dictBuilder/zdict.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zstd-0.7.4/dictBuilder/divsufsort.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zstd-0.7.4/decompress/zstd_decompress.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zstd-0.7.4/decompress/huf_decompress.o",
                    "./frontend/thirdparty/blosc/internal-complibs/zstd-0.7.4/decompress/zbuff_decompress.o",
                    "./frontend/thirdparty/blosc/internal-complibs/snappy-1.1.1/snappy-c.o",
                    "./frontend/thirdparty/blosc/internal-complibs/snappy-1.1.1/snappy.o",
                    "./frontend/thirdparty/blosc/internal-complibs/snappy-1.1.1/snappy-sinksource.o",
                    "./frontend/thirdparty/blosc/internal-complibs/snappy-1.1.1/snappy-stubs-internal.o",
                ],
                # library_dirs = [
                #         "/usr/lib64/openmpi/lib/",
                # ],
                # libraries=[
                #         "mpi",
                # ],
            ),
        ]),
        include_dirs = [
                # "/usr/include/openmpi-x86_64/",
                # "./thirdparty/blosc",
        ],
)

# python setup.py build_ext --inplace

