const std = @import("std");
const Allocator = std.mem.Allocator;
const testing = std.testing;

const BytePairError = error{
    NonAsciiInput,
    NullByteInput,
    UnrecognisedIndx,
    IndxTooBig,
    TruncatedSymbol,
    CorpusTooBig,
};

fn isAscii(c: u8) bool {
    return c & 0x80 == 0;
}

/// Variable length symbol returned by the read and write iterators
const Symbol = union(enum) {
    end: void,
    ascii: u7,
    indx: u15,

    pub fn init(raw: anytype) BytePairError!Symbol {
        if (@TypeOf(raw) == comptime_int) {
            if (raw < 0x80)
                return .{ .ascii = @intCast(raw) }
            else if (raw <= 0xffff)
                return .{ .indx = @truncate(raw) }
            else
                unreachable;
        } else if (@TypeOf(raw) == u8) {
            if (isAscii(raw))
                return .{ .ascii = @intCast(raw) }
            else
                return BytePairError.NonAsciiInput;
        }

        const bf: []const u8 = raw;
        const b0 = bf[0];

        if (isAscii(b0))
            return .{ .ascii = @intCast(b0) };

        if (bf.len < 2) {
            std.debug.print("truncated {}", .{b0});
            return BytePairError.TruncatedSymbol;
        }

        const b1 = bf[1];

        return .{ .indx = (@as(u15, b0) << 8) + b1 };
    }

    fn notEnd(self: Symbol) bool {
        return switch (self) {
            .end => false,
            else => true,
        };
    }

    fn order(self: Symbol) u16 {
        return switch (self) {
            .end => unreachable,
            .ascii => |i| @intCast(i),
            .indx => |i| 0x80 + @as(u16, @intCast(i)),
        };
    }

    fn eql(self: Symbol, b: Symbol) bool {
        return self.order() == b.order();
    }
};

/// A double buffered symbol buffer
const SymbolBuf = struct {
    side: u1 = 0,
    buf: []u8,
    len: usize,

    pub fn init(in: []const u8, buf: []u8) SymbolBuf {
        std.debug.assert(buf.len / 2 == in.len);

        @memcpy(buf[0..in.len], in);

        return SymbolBuf{
            .buf = buf,
            .len = in.len,
        };
    }

    fn getSide(self: SymbolBuf, side: u1) []u8 {
        const s = side * (self.buf.len / 2);

        return self.buf[s..][0..self.len];
    }

    fn getEncoded(self: SymbolBuf) []u8 {
        return self.getSide(self.side);
    }

    fn reader(self: SymbolBuf) SymbolReader {
        return SymbolReader{ .buf = self.getSide(self.side) };
    }

    fn writer(self: SymbolBuf) SymbolWriter {
        return SymbolWriter{ .buf = self.getSide(self.side ^ 1) };
    }

    fn flip(self: *SymbolBuf, written: SymbolWriter) void {
        self.len = written.i;
        self.side ^= 1;
    }
};

/// Iterate over symbols in a buffer
const SymbolReader = struct {
    i: usize = 0,
    buf: []const u8,

    fn next(self: *SymbolReader) BytePairError!Symbol {
        const ln = self.buf.len;

        if (self.i >= ln)
            return Symbol.end;

        const bf = self.buf[self.i..];
        const sym = try Symbol.init(bf);

        self.i += switch (sym) {
            .end => 0,
            .ascii => 1,
            .indx => 2,
        };

        return sym;
    }
};

/// Write symbols to a buffer
const SymbolWriter = struct {
    i: usize = 0,
    buf: []u8,

    fn next(self: *SymbolWriter, sym: Symbol) void {
        const i = self.i;
        const buf = self.buf;

        switch (sym) {
            .ascii => |c| {
                buf[i] = c;

                self.i += 1;
            },
            .indx => |indx| {
                buf[i] = 0x80 + @as(u8, @intCast(indx >> 8));
                buf[i + 1] = @truncate(indx);

                self.i += 2;
            },
            .end => {},
        }
    }

    fn rest(self: *SymbolWriter, reader: *SymbolReader) void {
        const left = reader.buf.len - reader.i;

        @memcpy(self.buf[self.i..][0..left], reader.buf[reader.i..]);

        reader.i += left;
        self.i += left;
    }
};

test "SymbolBuf ASCII" {
    const corpus = "abc DEF";
    var buf: [2 * corpus.len]u8 = undefined;

    var syms = SymbolBuf.init(corpus, &buf);

    for (0..2) |_| {
        var read = syms.reader();

        var s0 = try read.next();
        var s1 = try read.next();
        var i: usize = 0;

        try testing.expectEqual(@as(usize, 2), read.i);

        while (s1.notEnd()) : ({
            s0 = s1;
            s1 = try read.next();
            i += 1;
        }) {
            try testing.expectEqual(corpus[i], s0.ascii);
            try testing.expectEqual(corpus[i + 1], s1.ascii);
        }

        try testing.expectEqual(corpus.len - 1, i);

        read = syms.reader();
        var write = syms.writer();

        s0 = try read.next();

        while (s0.notEnd()) : (s0 = try read.next()) {
            write.next(s0);
        }

        syms.flip(write);
    }
}

test "SymbolBuf Index" {
    const corpus = [_]u8{ 'a', ' ', 0x80, 0x00, 0xff, 0xff };
    var buf: [2 * corpus.len]u8 = undefined;

    var syms = SymbolBuf.init(&corpus, &buf);

    for (0..2) |_| {
        var rdr = syms.reader();

        var s0 = try rdr.next();
        var s1 = try rdr.next();

        const Si = Symbol.init;
        try testing.expectEqual(try Si('a'), s0);
        try testing.expectEqual(try Si(' '), s1);

        s0 = s1;
        s1 = try rdr.next();

        try testing.expectEqual(try Si(' '), s0);
        try testing.expectEqual(try Si(0x8000), s1);

        s0 = s1;
        s1 = try rdr.next();

        try testing.expectEqual(try Si(0x8000), s0);
        try testing.expectEqual(try Si(0xffff), s1);

        try testing.expectEqual(Symbol.end, try rdr.next());

        rdr = syms.reader();
        var write = syms.writer();

        s0 = try rdr.next();

        while (s0.notEnd()) : (s0 = try rdr.next()) {
            write.next(s0);
        }

        syms.flip(write);
    }
}

const SymbolPairHashCtx = struct {
    const Self = @This();

    pub fn hash(self: Self, key: [2]Symbol) u32 {
        _ = self;
        const k = (@as(u32, key[0].order()) << 16) + key[1].order();

        return std.hash.uint32(k);
    }

    pub fn eql(self: Self, a: [2]Symbol, b: [2]Symbol, b_index: usize) bool {
        _ = self;
        _ = b_index;

        return a[0].eql(b[0]) and a[1].eql(b[1]);
    }
};

fn SymbolPairHashMap(comptime V: type) type {
    return std.ArrayHashMap([2]Symbol, V, SymbolPairHashCtx, false);
}

/// Byte pair transcoder
pub const BytePair = struct {
    const tableLen = std.math.maxInt(u15);

    alloc: Allocator,
    minFreq: u32,
    freq: SymbolPairHashMap(u32),
    enc: SymbolPairHashMap(Symbol) = undefined,
    dec: [tableLen][2]Symbol = undefined,
    maxIndx: u15 = 0,

    /// Create a transcoder for the given corpus
    /// If the corpus is very large then this is likely to be unbearably slow
    pub fn init(alloc: Allocator, minFreq: u32) !BytePair {
        var self = BytePair{
            .alloc = alloc,
            .minFreq = minFreq,
            .freq = SymbolPairHashMap(u32).init(alloc),
            .enc = SymbolPairHashMap(Symbol).init(alloc),
        };

        try self.enc.ensureTotalCapacity(tableLen);
        try self.freq.ensureTotalCapacity(1024);

        return self;
    }

    fn deinit(self: *BytePair) void {
        self.enc.clearAndFree();
        self.enc.deinit();
        self.maxIndx = 0;
        self.freq.clearAndFree();
        self.freq.deinit();
    }

    fn addCorpus(self: *BytePair, corpus: []const u8) !void {
        var buf: [2048]u8 = undefined;

        if (corpus.len > buf.len / 2)
            return BytePairError.CorpusTooBig;

        for (corpus) |c| {
            if (!isAscii(c))
                return BytePairError.NonAsciiInput;
        }

        var syms = SymbolBuf.init(corpus, buf[0 .. 2 * corpus.len]);
        _ = try self.encodeInplace(&syms);

        while (self.maxIndx < tableLen) {
            var maxFreq: u32 = 0;
            var mostFreq = [_]Symbol{Symbol{ .end = {} }} ** 2;

            var read = syms.reader();
            var s0 = try read.next();
            var s1 = try read.next();

            while (s1 != .end) : ({
                s0 = s1;
                s1 = try read.next();
            }) {
                const f = try self.freq.getOrPut(.{ s0, s1 });

                if (!f.found_existing)
                    f.value_ptr.* = 0;

                f.value_ptr.* += 1;
                const fv = f.value_ptr.*;

                if (maxFreq < fv) {
                    maxFreq = fv;
                    mostFreq[0] = s0;
                    mostFreq[1] = s1;
                }
            }

            if (maxFreq < self.minFreq)
                break;

            self.freq.clearRetainingCapacity();

            self.dec[self.maxIndx][0] = mostFreq[0];
            self.dec[self.maxIndx][1] = mostFreq[1];

            const substSym = Symbol{ .indx = self.maxIndx };
            const e = try self.enc.getOrPut(.{ mostFreq[0], mostFreq[1] });
            //std.debug.print("\nmf {} {} {} {}\n", .{ maxFreq, mostFreq[0], mostFreq[1], substSym });
            std.debug.assert(!e.found_existing);
            e.value_ptr.* = substSym;

            self.maxIndx += 1;

            read = syms.reader();
            s0 = try read.next();
            s1 = try read.next();

            var write = syms.writer();
            var subs: u32 = 0;

            while (true) {
                if (s0.eql(mostFreq[0]) and s1.eql(mostFreq[1])) {
                    //std.debug.print("subst {} {} -> {}\n", .{ s0, s1, substSym });

                    write.next(substSym);

                    subs += 1;
                    if (subs >= maxFreq)
                        break;

                    s0 = try read.next();
                    s1 = try read.next();
                } else {
                    //std.debug.print("keep {}\n", .{s0});

                    write.next(s0);

                    s0 = s1;
                    s1 = try read.next();
                }

                if (s1 == .end) {
                    write.next(s0);
                    break;
                }
            }

            write.rest(&read);
            syms.flip(write);
        }
    }

    fn encodeInplace(self: BytePair, syms: *SymbolBuf) ![]u8 {
        while (true) {
            var any: bool = false;
            var read = syms.reader();
            var s0 = try read.next();
            var s1 = try read.next();
            var r: u15 = 0x7fff;
            var last: usize = 0;

            while (s1 != .end) : ({
                s0 = s1;
                s1 = try read.next();
            }) {
                if (self.enc.get(.{ s0, s1 })) |sym| {
                    switch (sym) {
                        .indx => |i| {
                            any = true;

                            if (i < r) {
                                r = i;
                                last = read.i;
                            }
                        },
                        else => unreachable,
                    }
                }
            }

            if (!any)
                break;

            const rs = Symbol{ .indx = r };
            const rs0 = self.dec[r][0];
            const rs1 = self.dec[r][1];

            var write = syms.writer();
            read = syms.reader();

            s0 = try read.next();
            s1 = try read.next();

            while (true) {
                if (rs0.eql(s0) and rs1.eql(s1)) {
                    //std.debug.print("subst {} {} -> {}\n", .{ s0, s1, rs });
                    write.next(rs);

                    if (read.i >= last)
                        break;

                    s0 = try read.next();
                    s1 = try read.next();
                } else {
                    //std.debug.print("keep {}\n", .{s0});

                    write.next(s0);

                    s0 = s1;
                    s1 = try read.next();
                }

                if (s1 == .end) {
                    write.next(s0);
                    break;
                }
            }

            write.rest(&read);

            syms.flip(write);
        }

        return syms.getSide(syms.side);
    }

    /// Encode a slice of bytes allocating the output
    fn encodeAlloc(self: BytePair, in: []const u8) ![]u8 {
        const bsize = 1024;
        var buf: [2 * bsize]u8 = undefined;
        const out = try self.alloc.alloc(u8, in.len);

        var i: usize = 0;
        var o: usize = 0;
        while (i + bsize <= in.len) : (i += bsize) {
            var syms = SymbolBuf.init(in[i..][0..bsize], &buf);
            const encoded = try self.encodeInplace(&syms);
            @memcpy(out[o..][0..syms.len], encoded);
            o += syms.len;
        }
        if (i < in.len) {
            var syms = SymbolBuf.init(in[i..], buf[0 .. 2 * (in.len - i)]);
            const encoded = try self.encodeInplace(&syms);
            @memcpy(out[o..][0..syms.len], encoded);
            o += syms.len;
        }

        return try self.alloc.realloc(out, o);
    }

    /// Decode a slice of bytes allocating the output
    fn decodeAlloc(self: BytePair, in: []const u8) ![]u8 {
        var stack = try std.ArrayList(Symbol).initCapacity(
            self.alloc,
            @max(1, std.math.log2(self.maxIndx)),
        );
        defer stack.deinit();
        var out = try std.ArrayList(u8).initCapacity(self.alloc, 2 * in.len);
        defer out.deinit();
        var read = SymbolReader{ .buf = in };
        var s = try read.next();

        while (s != .end) {
            switch (s) {
                .end => unreachable,
                .ascii => |c| {
                    try out.append(c);

                    //std.debug.print("add {}\n", .{c});

                    s = if (stack.popOrNull()) |sym|
                        sym
                    else
                        try read.next();
                },
                .indx => |i| {
                    if (i > self.maxIndx)
                        return BytePairError.UnrecognisedIndx;

                    s = self.dec[i][0];

                    try stack.append(self.dec[i][1]);

                    //std.debug.print("i {} d {} p {}\n", .{ i, s, self.dec[i][1] });
                },
            }
        }

        return try out.toOwnedSlice();
    }
};

test "BytePair alpha" {
    const alpha = "abcdefghijklmnopqrstuwxyz1";

    var bp = try BytePair.init(testing.allocator, 1);
    defer bp.deinit();

    try bp.addCorpus(alpha);

    const encoded = try bp.encodeAlloc(alpha);
    defer testing.allocator.free(encoded);
    try testing.expectEqual(@as(usize, 2), encoded.len);
    try testing.expectEqual(bp.maxIndx - 1, (@as(u15, @intCast(encoded[0])) << 8) + encoded[1]);

    const decoded = try bp.decodeAlloc(encoded);
    defer testing.allocator.free(decoded);

    try testing.expectEqualStrings(alpha, decoded);
}

test "BytePair sentences" {
    const alpha = "The quick brown fox jumped over the lazy white dog. The slow badger made a run for the hotel bar, but in the event was too late and went to be completely sober.";

    var bp = try BytePair.init(testing.allocator, 1);
    defer bp.deinit();

    try bp.addCorpus(alpha);

    const encoded = try bp.encodeAlloc(alpha);
    defer testing.allocator.free(encoded);
    try testing.expectEqual(@as(usize, 2), encoded.len);
    try testing.expectEqual(bp.maxIndx - 1, (@as(u15, @intCast(encoded[0])) << 8) + encoded[1]);

    const decoded = try bp.decodeAlloc(encoded);
    defer testing.allocator.free(decoded);
    try testing.expectEqualStrings(alpha, decoded);

    const alpha2 = "An unrelated sentence which maybe shares some sequences. We'll see how well it can be compressed.";

    const encoded2 = try bp.encodeAlloc(alpha2);
    defer testing.allocator.free(encoded2);

    std.debug.print("\ncompression {} / {} = {}\n", .{ encoded2.len, alpha2.len, @as(f64, @floatFromInt(encoded2.len)) / @as(f64, @floatFromInt(alpha2.len)) });

    const decoded2 = try bp.decodeAlloc(encoded2);
    defer testing.allocator.free(decoded2);
    try testing.expectEqualStrings(alpha2, decoded2);
}

test "BytePair ASCII" {
    const ascii = init: {
        var arr: [2 * 128 * 128]u8 = undefined;

        var k: usize = 0;

        for (0..128) |i| {
            for (0..128) |j| {
                arr[k] = @intCast(i);
                arr[k + 1] = @intCast(j);
                k += 2;
            }
        }

        if (k != arr.len) unreachable;

        break :init arr;
    };

    var bp = try BytePair.init(testing.allocator, 2);
    defer bp.deinit();

    var r: usize = 0;
    while (r < ascii.len) : (r += 1024) {
        const block = ascii[r .. r + 1024];

        try bp.addCorpus(block);
    }
    std.debug.print("Dictionary size {}", .{bp.maxIndx});

    const encoded = try bp.encodeAlloc(&ascii);
    defer testing.allocator.free(encoded);
    std.debug.print("compression {} / {} = {}", .{ encoded.len, ascii.len, @as(f64, @floatFromInt(encoded.len)) / @as(f64, @floatFromInt(ascii.len)) });
    const decoded = try bp.decodeAlloc(encoded);
    defer testing.allocator.free(decoded);

    try testing.expectEqualStrings(&ascii, decoded);
}
