import sys
import gio

name = sys.argv[1]
gio.gio_inspect(name)

x = gio.gio_read(name, "x")
print x

