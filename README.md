# GenericIO

Forked from https://trac.alcf.anl.gov/projects/genericio

The readme is my attempt to understand the code. This is complicated by the fact that I don't know any C++ :)

## Questions
* What is the `MismatchBehavior` in the header reading funcs.
* What does the `EffRank` variable mean.
* What are the variable flags for either the extra space or the IsPhysCoordX
    * Why do we need the extra space on read? Why is it always just 8.
* Can you append a column?
* Using structured arrays makes a nice API but is horrible for performance because we keep needing to switch from column contiguous (gio) to row contiguous (numpy).
* When reading with a different number of ranks to what was written, how should we partition ranks? There must be a way to easily get all nearby ranks.
  * I suppose that this depends on how they were partitioned before. On read you should probably partition so that the ranks that wrote to the same file are read by the same rank?
* Why is performance so much worse when reading the coreproperties vs the test data here.
  * coreproperties was written with many more ranks (256 vs 8)
  * coreproperties has many more variables (in the tests in both cases we only read 2)

## Partition

We sometimes split COMM_WORLD into smaller groups. The partition arg determines which of these groups you belong in. By default there is just the global partition.

The reason to do this is IO performance. Ideally you want to stripe your data over all of the storage devices on the computer. E.g. if you have 10 disks and a 100MB file, write 10MB to each disk. There might be a limit to the the number of drives you can stripe a single file over (or just a default). In that case you can break up your file into smaller files each of which will be striped.

The example in the gio docs says: You want the number of files * the number of stripes per file to be roughly equal to the number of drives in the storage system. So for titan with ~1k OSTs (object storage targets) and 4 stripes per file you want to be writing ~250 files. See the [titan file system docs](https://www.olcf.ornl.gov/for-users/system-user-guides/titan/file-systems/) which describes Lustre (the file system) quite nicely.

## BigEndian

Lots of references to this in the code. I think it is to handle the case where data is written by a little endian machine and read by a big endian one (or vice versa). Naively, that wouldn't work as the little endian machine would write the integer 1 as `01 00 00 00` which the big endian machine would read as 2^24.

To fix this, we byteswap the data if the endianness isn't the same in the data and on our machine (I think).

## Blosc

This uses [blosc](http://blosc.org/pages/blosc-in-depth/) which reads and writes compressed data. This not only results in smaller disk/ram usage but is also faster! The cost of decompressing is smaller than the cost of loading the larger data set. see [a paper](http://blosc.org/pages/blosc-in-depth/) for some more details.

To understand Gio, all you need to know is that blosc is lowest layer in reading/writing and it does it fast.

## GenericFileIO

This isn't that important but I didn't know that when I first looked at this...

Base class `GenericFileIO` with two subclasses `GenericFileIO_MPI` and `GenericFileIO_POSIX`. There is also a `GenericFileIO_MPICollective` which is based off `GenericFileIO_MPI`. These are simple readers/writers. But that are coordinated by the `GenericIO` class.

#### GenericFileIO_POSIX

```
class GenericFileIO_POSIX : public GenericFileIO {
public:
  GenericFileIO_POSIX() : FH(-1) {}
  ~GenericFileIO_POSIX();

public:
  void open(const std::string &FN, bool ForReading = false);
  void setSize(size_t sz);
  void read(void *buf, size_t count, off_t offset, const std::string &D);
  void write(const void *buf, size_t count, off_t offset, const std::string &D);

protected:
  int FH;
};

class GenericFileIO {
// Ignored
protected:
  std::string FileName;
};
```

* FH: File handle/file descriptor

* Constructor: Initialized FH to -1 (uninitilized)
* Destructor: Closes FH if if it is not uninitialized

* open: Sets `FH` to the file descriptor given by opening the file. If `ForReading`, open read only. Else open readwrite. Sets filename in the base class.
* setSize: Truncate file to given length. File must be opened for writing else this fails.
* read: Read `count` bytes into `buf` starting at byte `offset`. `D` is added to the error string on error.
* write: Same as read, just in reverse.

#### GenericFileIO_MPI

Appears to be pretty much like the POSIX version except using MPI calls (e.g. `MPI_file_open` vs `open`)? Open, close, and set size are collective. Write and read aren't.

#### GenericFileIO_MPICollective

Subclass of GenericFileIO_MPI where read and write are collective.
