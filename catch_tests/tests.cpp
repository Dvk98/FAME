#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "ShiftAnd.h"
#include "ReadQueue.h"
#include "BitFunctions.h"

std::array<uint8_t, 16> Setup() {
    std::array<uint8_t, 16> lmap;
    lmap['A'%16] = 0;
    lmap['C'%16] = 1;
    lmap['G'%16] = 2;
    lmap['T'%16] = 3;
    return lmap;
}


TEST_CASE("get Highest Index for uint_64", "[getHighestIdx64]") {
    uint64_t n = 0;
    REQUIRE(BitFun::getHighestIdx64(n) == 0);
    n = 64;
    REQUIRE(BitFun::getHighestIdx64(n) == 7);
}

TEST_CASE("Local ShiftAnd_Test", "[localShiftAnd]") {
    std::array<uint8_t, 16> lmap = Setup();
    std::string read   = "AAAAAAAA";
    std::string window = "AAAAAAAT";
    std::vector<char> t (window.begin(), window.end());

    ShiftAnd<0> sa0(read, lmap);
    std::vector<uint64_t> matchings0;
    std::vector<uint8_t> errors0;
    std::vector<uint8_t> length0;
    sa0.querySeqLocal(t.begin(), t.end(), matchings0, errors0, length0, 4);

    REQUIRE(length0.size() == 4);
    REQUIRE(length0[0] == 4);
    REQUIRE(length0[1] == 5);
    REQUIRE(length0[2] == 6);
    REQUIRE(length0[3] == 7);

    REQUIRE(matchings0.size() == 4);
    REQUIRE(matchings0[0] - (length0[0] - 1) == 0);
    REQUIRE(matchings0[1] - (length0[1] - 1) == 0);
    REQUIRE(matchings0[2] - (length0[2] - 1) == 0);
    REQUIRE(matchings0[3] - (length0[3] - 1) == 0);

    REQUIRE(errors0.size() == 4);
    REQUIRE(errors0[0] == 0);
    REQUIRE(errors0[1] == 0);
    REQUIRE(errors0[2] == 0);
    REQUIRE(errors0[3] == 0);
}

TEST_CASE("Local ShiftAnd_Test with splitted matches", "[localShiftAnd]") {
    std::array<uint8_t, 16> lmap = Setup();
    std::string read   = "AAAAAAAA";
    std::string window = "AAAAAAATAAAAAA";
    std::vector<char> t (window.begin(), window.end());

    ShiftAnd<0> sa0(read, lmap);
    std::vector<uint64_t> matchings0;
    std::vector<uint8_t> errors0;
    std::vector<uint8_t> length0;
    sa0.querySeqLocal(t.begin(), t.end(), matchings0, errors0, length0, 6);

    REQUIRE(length0.size() == 3);
    REQUIRE(length0[0] == 6);
    REQUIRE(length0[1] == 7);
    REQUIRE(length0[2] == 6);

    REQUIRE(matchings0.size() == 3);
    REQUIRE(matchings0[0] - (length0[0] - 1) == 0);
    REQUIRE(matchings0[1] - (length0[1] - 1) == 0);
    REQUIRE(matchings0[2] - (length0[2] - 1) == 8);

    REQUIRE(errors0.size() == 3);
    REQUIRE(errors0[0] == 0);
    REQUIRE(errors0[1] == 0);
    REQUIRE(errors0[2] == 0);
}

TEST_CASE("Local ShiftAnd_Test_Reverse", "[localShiftAnd]") {
    std::array<uint8_t, 16> lmap = Setup();
    std::string read   = "TTTTTTTT";
    std::string window = "TAAAAAAA"; //AAAAAAAT -> TTTTTTTA -> ATTTTTTT
    std::vector<char> t (window.begin(), window.end());

    ShiftAnd<0> sa0(read, lmap);
    std::vector<uint64_t> matchings0;
    std::vector<uint8_t> errors0;
    std::vector<uint8_t> length0;
    sa0.queryRevSeqLocal(t.end() - 1, t.begin() - 1, matchings0, errors0, length0, 4);

    REQUIRE(length0.size() == 4);
    REQUIRE(length0[0] == 4);
    REQUIRE(length0[1] == 5);
    REQUIRE(length0[2] == 6);
    REQUIRE(length0[3] == 7);

    REQUIRE(matchings0.size() == 4);
    REQUIRE(matchings0[0] == 7);
    REQUIRE(matchings0[1] == 7);
    REQUIRE(matchings0[2] == 7);
    REQUIRE(matchings0[3] == 7);

    REQUIRE(errors0.size() == 4);
    REQUIRE(errors0[0] == 0);
    REQUIRE(errors0[1] == 0);
    REQUIRE(errors0[2] == 0);
    REQUIRE(errors0[3] == 0);
}


TEST_CASE("saQuerySeedSetRefLocal_Test", "[saQuerySeedSetRefLocal]") {
    std::array<uint8_t, 16> lmap = Setup();

    RefGenome ref("/home/dvk/FAME/catch_tests/RefGenomeTests/indexSmall");
    ReadQueue rq ("test.fastq", ref, false, false);
    rq.localAlign = true;

    std::string read = "CCACCTGATAATAAATATGTGTCAGAGAGGACTCGATCCGTGTGACCTTGGGTTACAGGGGTACTGGACCGTGAAAAGCCTTTGACAAAATCCTTTTGCC"; // 6th row + 9th row
    ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS> sa(read, lmap);
    MATCH::match match = 0;
    uint8_t length = 0;
    uint16_t qThreshold = 5;
    uint8_t minlength = 30;
    rq.getSeedRefs(read, read.size(), qThreshold);
    rq.saQuerySeedSetRefLocal(sa, match, length, qThreshold, minlength);
    std::cout << length << std::endl;
}


//TEST_CASE("matchLocalAlign_Test", "[matchLocalAlign]") {
//
//}

