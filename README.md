# GenericIO

Forked from https://trac.alcf.anl.gov/projects/genericio

The readme is my attempt to understand the code. This is complicated by the fact that I don't know any C++ :)

## Interface

Base class `GenericFileIO` with two subclasses `GenericFileIO_MPI` and `GenericFileIO_POSIX`. There is also a `GenericFileIO_MPICollective` which is based off `GenericFileIO_MPI`. I think these are simple readers/writers. But that are coordinated by the `GenericIO` class.

What does the GenericIO class do...

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
