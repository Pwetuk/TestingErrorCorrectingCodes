#!/usr/bin/env bash

# This will sum up the sizes of *all* sections in your ELF
# that are marked ALLOC+WRITE (i.e. actually loaded into RAM at startup).
readelf -W -S my_program \
  | awk '
    $0 ~ /\[ *[0-9]+\] / {
      name = $2;      # section name
      flags = $4;     # flags string, e.g. "WA" or "AX"
      size  = strtonum("0x"$6);
      # If ALLOC and WRITE but not NOBITS-only (we want PROGBITS too)
      if (flags ~ /A/ && flags ~ /W/) {
        alloc_write += size;
      }
      # If ALLOC but *not* WRITE (that's .rodata etc.), skip here
      # If NOBITS (BSS), we'll count that separately below.
    }
    END {
      # Now also grab BSS explicitly (ALLOC but *no* LOADed data):
      # BSS is ALLOC + WRITE but NOBITS, but readelf shows size in the same table.
      # Actually the above test already picks up BSS (it’s ALLOC+WRITE), so
      # alloc_write includes both DATA+NOBITS.
      printf("static RAM (data + bss) = %d bytes\n", alloc_write);
    }'
