const std = @import("std");
const testing = std.testing;
const enc = @import("enc.zig");

test {
    testing.refAllDecls(enc);
}
