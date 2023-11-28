const std = @import("std");
const Allocator = std.mem.Allocator;
const testing = std.testing;

const BytePairError = error{
    NonAsciiInput,
    IndxTooBig,
    TruncatedSymbol,
    CorpusTooBig,
};

fn isAscii(c: u8) bool {
    return c & 0x80 == 0;
}

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

        if (bf.len < 2)
            return BytePairError.TruncatedSymbol;

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
            .indx => |i| 0x7f + @intCast(i),
        };
    }

    fn eql(self: Symbol, b: Symbol) bool {
        return self.order() == b.order();
    }
};

const SymbolBuf = struct {
    side: u1 = 0,
    buf: []u8,

    pub fn init(corpus: []const u8, buf: []u8) SymbolBuf {
        std.debug.assert(buf.len / 2 == corpus.len);

        @memcpy(buf[0..corpus.len], corpus);

        return SymbolBuf{
            .buf = buf,
        };
    }

    fn len(self: SymbolBuf) usize {
        return self.buf.len / 2;
    }

    fn getSide(self: SymbolBuf, side: u1) []u8 {
        const ln = self.len();
        const s = side * ln;

        return self.buf[s..][0..ln];
    }

    fn reader(self: SymbolBuf) SymbolReader {
        return SymbolReader{ .buf = self.getSide(self.side) };
    }

    fn writer(self: SymbolBuf) SymbolWriter {
        return SymbolWriter{ .buf = self.getSide(self.side ^ 1) };
    }

    fn flip(self: *SymbolBuf) void {
        self.side ^= 1;
    }
};

const SymbolReader = struct {
    i: usize = 0,
    buf: []u8,

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
                buf[i] = @intCast(indx >> 7);
                buf[i + 1] = @truncate(indx);

                self.i += 2;
            },
            .end => unreachable,
        }
    }
};

test "SymbolBuf ASCII" {
    const corpus = "abc DEF";
    var buf: [2 * corpus.len]u8 = undefined;

    var syms = SymbolBuf.init(corpus, &buf);

    for (0..1) |_| {
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

        read = syms.reader();
        var write = syms.writer();

        s0 = try read.next();

        while (s0.notEnd()) : (s0 = try read.next()) {
            write.next(s0);
        }

        syms.flip();
    }
}

test "SymbolBuf Index" {
    const corpus = [_]u8{ 'a', ' ', 0x80, 0x01, 0xff, 0xff };
    var buf: [2 * corpus.len]u8 = undefined;

    var syms = SymbolBuf.init(&corpus, &buf);
    var rdr = syms.reader();

    var s0 = try rdr.next();
    var s1 = try rdr.next();

    const Si = Symbol.init;
    try testing.expectEqual(try Si('a'), s0);
    try testing.expectEqual(try Si(' '), s1);

    s0 = s1;
    s1 = try rdr.next();

    try testing.expectEqual(try Si(' '), s0);
    try testing.expectEqual(try Si(0x8001), s1);

    s0 = s1;
    s1 = try rdr.next();

    try testing.expectEqual(try Si(0x8001), s0);
    try testing.expectEqual(try Si(0xffff), s1);

    try testing.expectEqual(Symbol.end, try rdr.next());
}

const SymbolPairHashCtx = struct {
    const Self = @This();

    pub fn hash(_: Self, key: [2]Symbol) u32 {
        const k = (@as(u32, key[0].order()) << 16) + [1]key.order();

        return std.hash.uint32(k);
    }

    pub fn eql(_: Self, a: [2]Symbol, b: [2]Symbol) bool {
        return a[0].eql(b[0]) and a[1].eql(b[1]);
    }
};

fn SymbolPairHashMap(comptime V: type) type {
    return std.HashMap([2]Symbol, V, SymbolPairHashCtx, 80);
}

pub const BytePair = struct {
    const tableLen = std.math.maxInt(u15);

    enc: SymbolPairHashMap(Symbol) = undefined,
    dec: [tableLen][2]Symbol = undefined,
    maxIndx: u15 = 0,

    pub fn init(alloc: Allocator, corpus: []const u8) !BytePair {
        for (corpus) |c| {
            if (!isAscii(c))
                return BytePairError.NonAsciiInput;
        }

        if (corpus.len > std.math.maxInt(u32))
            return BytePairError.CorpusTooBig;

        const self = BytePair{
            .enc = SymbolPairHashMap(Symbol).init(alloc),
        };

        try self.enc.ensureTotalCapacity(tableLen);

        const bufs = try alloc.alloc(u8, 2 * corpus.len);
        defer alloc.free(bufs);
        var syms = SymbolBuf.init(corpus, bufs);

        var freq = SymbolPairHashMap(u32).init(alloc);
        defer freq.deinit();
        try freq.ensureTotalCapacity(@min(corpus.len, tableLen));

        while (self.maxIndx < tableLen) {
            var maxFreq: u32 = 0;
            var mostFreq = .{Symbol{.end}} ** 2;

            var read = syms.reader();
            var s0 = read.next();
            var s1 = read.next();

            while (s1 != .end) : ({
                s0 = s1;
                s1 = read.next();
            }) {
                const f = try freq.getOrPut(.{ s0, s1 });

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

            if (maxFreq < 1)
                break;

            freq.clearRetainingCapacity();

            self.dec[self.maxIndx][0] = s0;
            self.dec[self.maxIndx][1] = s1;

            const substSym = Symbol{ .indx = self.maxIndx };
            const e = try self.enc.getOrPut(.{ mostFreq[0], mostFreq[1] });
            std.debug.assert(!e.found_existing);
            e.value_ptr.* = substSym;

            self.maxIndx += 1;

            read = syms.reader();
            s0 = read.next();
            s1 = read.next();

            var write = syms.writer();

            while (s1 != .end) : ({
                s0 = s1;
                s1 = read.next();
            }) {
                if (s0.eql(mostFreq[0]) and s1.eql(mostFreq[1])) {
                    write.next(substSym);
                } else {
                    write.next(s0);
                }
            }

            write.next(s0);
            write.next(Symbol{.end});
            syms.flip();
        }

        return self;
    }
};

test "BytePair" {}
