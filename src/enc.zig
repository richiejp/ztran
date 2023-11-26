const std = @import("std");
const testing = std.testing;

const BytePairError = error{
    NonAsciiInput,
    IndxTooBig,
    TruncatedSymbol,
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

pub const BytePair = struct {
    const dicLen = std.math.maxInt(u16) >> 1;

    dic: [dicLen]u16 = undefined,

    pub fn init(corpus: []const u8) !BytePair {
        const self = BytePair{};
        _ = self;

        for (corpus) |c| {
            if (!isAscii(c))
                return BytePairError.NonAsciiInput;
        }

        const bufs: [2 * corpus.len]u8 = undefined;
        var syms = SymbolBuf.init(&corpus, &bufs);
        var reader = syms.reader();
        const freqMat = [_]([]u16){.{0} ** dicLen} ** dicLen;
        _ = freqMat;

        var s0 = rdr.next();
        var s1 = rdr.next();

        while (s1 != .end) : ({
            s0 = rdr.next();
            s1 = rdr.next();
        }) {}
    }
};
