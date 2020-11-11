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
    REQUIRE(matchings0[0] == 0);
    REQUIRE(matchings0[1] == 0);
    REQUIRE(matchings0[2] == 0);
    REQUIRE(matchings0[3] == 0);

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
    REQUIRE(matchings0[0] == 0);
    REQUIRE(matchings0[1] == 0);
    REQUIRE(matchings0[2] == 8);

    REQUIRE(errors0.size() == 3);
    REQUIRE(errors0[0] == 0);
    REQUIRE(errors0[1] == 0);
    REQUIRE(errors0[2] == 0);
}

TEST_CASE("Local ShiftAnd_Test_Reverse", "[localShiftAnd]") {
    std::array<uint8_t, 16> lmap = Setup();
    std::string read   = "TTTTTTTT";
    std::string window = "AAAAAAAT";
    std::vector<char> t (window.begin(), window.end());

    ShiftAnd<0> sa0(read, lmap);
    std::vector<uint64_t> matchings0;
    std::vector<uint8_t> errors0;
    std::vector<uint8_t> length0;
    sa0.queryRevSeqLocal(t.begin(), t.end(), matchings0, errors0, length0, 4);

    REQUIRE(length0.size() == 4);
    REQUIRE(length0[0] == 4);
    REQUIRE(length0[1] == 5);
    REQUIRE(length0[2] == 6);
    REQUIRE(length0[3] == 7);

    REQUIRE(matchings0.size() == 4);
    REQUIRE(matchings0[0] == 0);
    REQUIRE(matchings0[1] == 0);
    REQUIRE(matchings0[2] == 0);
    REQUIRE(matchings0[3] == 0);

    REQUIRE(errors0.size() == 4);
    REQUIRE(errors0[0] == 0);
    REQUIRE(errors0[1] == 0);
    REQUIRE(errors0[2] == 0);
    REQUIRE(errors0[3] == 0);
}


TEST_CASE("saQuerySeedSetRefLocal_Test", "[saQuerySeedSetRefLocal]") {
    std::array<uint8_t, 16> lmap = Setup();
    // set up sequence container
    std::string seq       = "ATGTTGCCTAATTTCACTATTCAGGGTTATACGCCTGGAATATTCTAGGATTCCTAGTCAATTTAT";
    // sequence with reduced alphabet
    std::string redSeq    = "ATGTTGTTTAATTTTATTATTTAGGGTTATATGTTTGGAATATTTTAGGATTTTTAGTTAATTTAT";
    // reverse sequence
    std::string revSeq    = "ATAAATTGACTAGGAATCCTAGAATATTCCAGGCGTATAACCCTGAATAGTGAAATTAGGCAACAT";
    std::string redRevSeq = "ATAAATTGATTAGGAATTTTAGAATATTTTAGGTGTATAATTTTGAATAGTGAAATTAGGTAATAT";
    std::vector<char> seqV (seq.begin(), seq.end());
    std::vector<std::vector<char> > genSeq;
    genSeq.push_back(seqV);

    // set up CpG container
    std::vector<struct CpG> cpgTab;
    std::vector<struct CpG> cpgStart;
    cpgStart.push_back({0, 3});


    std::unordered_map<uint8_t, std::string> chrMap;
    std::cout << "Test1" << std::endl;
    RefGenome ref (std::move(cpgTab), std::move(cpgStart), genSeq);
    //RefGenome ref (genSeq);
    std::cout << "Test2" << std::endl;
    ReadQueue rq ("test.fa", ref, false, false);
    rq.localAlign = true;

    std::string read = "ATGTTGCCTAATTTCACTATTTATTCTAGGATTCCTAGTCAATTTAT";
    ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS> sa(read, lmap);
    MATCH::match match;
    uint8_t length;
    uint16_t qThreshold;
    rq.saQuerySeedSetRefLocal(sa, match, length, qThreshold, 10);
    std::cout << length << std::endl;
}


//TEST_CASE("matchLocalAlign_Test", "[matchLocalAlign]") {
//
//}

