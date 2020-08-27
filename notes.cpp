For the match it would be possible to create another Match construct or just have a uint8_t length getting passed along with the match, which records the length of the match (e.g when the match is ReadLen long then 0)
For implementation i have 2 possibilites in mind:
First one is to call the method saQuerySeedSetRefLocal two times when no match is found. 
The First time for the first partial match of the sequence and the second one for the second partial match of the sequence. 
If both are found then we have a valid match. This would look like this:

//no match found at all
        } else {
            ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS> saFwdFirst(r.seq, lmap); 
			// Why is it necessary to create an automata for the reverse strang aswell? Because we are already processing the sequence against the forward & reverse reference genome
            int firstLength = 0;
            succQueryFwd = saQuerySeedSetRefLocal(saFwdFirst, matchFwdFirst, firstLength);

            if(succQueryFwd) {
                ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS> saFwdSecond(r.seq.substr(firstLength - 1), lmap);
                int secondLength = 0;
                succQueryFwd = saQuerySeedSetRefLocal(saFwdSecond, matchFwdSecond, secondLength);
            }

            if(succQueryFwd == 1) {
                ++succMatchT;
                // Compare with Rev Match and ...
                // Then computeMethLvl
            }
            else if (succQueryFwd == -1 || succQueryRev == -1)
            {
                r.isInvalid = true;
                ++nonUniqueMatchT;

            } else {
                r.isInvalid = true;
                ++unSuccMatchT;
            }
        }


saQuerySeedSetRefLocal(ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa, Match& match1, Match& match2, int length) {

    std::array<uint8_t, MyConst::ADDMIS + MyConst::MISCOUNT + 1> multiMatch;
    std::array<uint8_t, MyConst::ADDMIS + MyConst::MISCOUNT + 1> lengthMatch;
	multiMatch.fill(0);
    lengthMatch.fill(0);
    std::array<MATCH::match, MyConst::ADDMIS + MyConst::MISCOUNT + 1> uniqueMatches;

    for(const auto& m : fwdMetaIDs_t)
    {
        ...
        std::vector<uint64_t> matchings; // position in mcpg
		std::vector<uint8_t> errors; // error count
        std::vector<uint8_t> length; // length of matched sequence
        sa.querySeqLocal(startIt, endIt, matchings, errors, length);

        size_t i = 0;

        if (matchings.size() > 0)
		{
			// compare chromosome and offset
			if (matchings[0] + ref.metaWindows[m.first].startPos == prevOff && ref.metaWindows[m.first].chrom == prevChr)
			{
				++i;
			}
		}
		// go through matching and see if we had such a match (with that many errors) before - if so,
		// return to caller reporting no match
		for (; i < matchings.size(); ++i)
		{

			// check if we had a match with that many errors before
			if (multiMatch[errors[i]] && lengthMatch[errors[i]] == length[i]])
			{

				MATCH::match& match_2 = uniqueMatches[errors[i]];
				const bool isStart = MATCH::isStart(match_2);
				const bool isFwd = MATCH::isFwd(match_2);
				// check if same k-mer (borders of meta CpGs)
				if (isFwd && !isStart && ref.metaWindows[MATCH::getMetaID(match_2)].startPos + MATCH::getOffset(match_2) == ref.metaWindows[m.first].startPos + matchings[i])
				{
					continue;

				} else {
					if (!errors[i])
					{
						return -1;
					}
					multiMatch[errors[i]] = 2;
				}


			} else {
				uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i], errors[i], 1, 0, m.first);
                lengthMatch[errors[i]] = length[i]
				multiMatch[errors[i]] = 1;
			}
		}
		if (matchings.size() > 0)
		{

			prevChr = ref.metaWindows[m.first].chrom;
			prevOff = ref.metaWindows[m.first].startPos + matchings[matchings.size() - 1];

		} else {

			prevChr = 0;
			prevOff = 0xffffffffffffffffULL;
		}
	}
	prevChr = 0;
	prevOff = 0xffffffffffffffffULL;
    }

    for(const auto& m : revMetaIDs_t)
    {
        ...
        std::vector<uint64_t> matchings; // position in mcpg
		std::vector<uint8_t> errors; // error count
        std::vector<uint8_t> length; // length of matched sequence
        sa.queryRevSeqLocal(startIt, endIt, matchings, errors, length);

        size_t i = 0;

        // Checks are the same as for forward strang
    }


    // Same as original
	for (size_t i = 0; i < multiMatch.size(); ++i)
	{
		if (multiMatch[i] == 0)
		{
			continue;
		}
		mat = uniqueMatches[i];
		if (multiMatch[i] > 1)
		{
			return -1;
		} else {

			return 1;
		}
	}
	return 0;
}

//Second thing would be to just do the work in saQuerySeedSetRefLocal(), but then it would be needed to bitshift the Automata to cut the sequence right.