#include "ParamPack.h"
#include "ngtypes/basetypes.h"
#include "ngtypes/tree.h"


namespace PARAMPACK{
  //  void ComposeTracedCurves(IDataPointer m_Source)
  //  {
  //      NG_CREATE_DYNAMIC_CAST(NeuronPopulation, curTree, m_Source);
        //if (curTree->GetTree().size() < 2lu) return;//First level block, such as neuron dendrites will not be composed.
  //      /*reconstruct tree structure in current block*/
  //      for (size_t i = 0; i < curTree->GetTree().back().size(); ++i) {
  //          double minDist = 10000.0, tmpDist;
  //          Vec5d *pt;
  //          VectorVec5d &curCurve = curTree->GetTree().back()[i];
  //          Vec5d &curHeadNode = curCurve.front();
  //          for (size_t j = 0; j < curTree->GetTree().back().size(); ++j) {
  //              if (&curCurve == &curTree->GetTree().back()[j]) continue;
  //              for (size_t ij = 0; ij < curTree->GetTree().back()[j].size(); ++ij) {
  //                  tmpDist = (curTree->GetTree().back()[j][ij].block(0, 0, 3, 1) - curHeadNode.block(0, 0, 3, 1)).norm();
  //                  if (tmpDist < minDist) {
  //                      minDist = tmpDist;
  //                      pt = &curTree->GetTree().back()[j][ij];
  //                  }
  //              }
  //          }
  //          if (minDist < 5.0 && minDist > 0.99) {
  //              VectorVec5d tmpVectorVec5d; tmpVectorVec5d.reserve(curCurve.size() + 1);
  //              tmpVectorVec5d.push_back(*pt);
  //              std::copy(curCurve.begin(), curCurve.end(), std::back_inserter(tmpVectorVec5d));
  //              curCurve.swap(tmpVectorVec5d);
  //          }
  //      }
  //      /*connect curves in current block to previous curves*/
  //      int blockInd = -1, curveInd = -1, headInd = -1;
  //      std::vector<int > blockIndList;
  //      std::vector<int> curveIndList;
  //      for (size_t i = 0; i < curTree->GetTree().back().size(); ++i) {
  //          VectorVec5d &curCurve = curTree->GetTree().back()[i];
  //          if (curCurve.empty()) continue;
  //          PARAMPACK::FindNearestCurvesEnd(m_Source, curCurve, blockInd, curveInd);
  //          if (blockInd == -1) {
  //              printf("There are curves individual.\n");
  //              continue;
  //          }
  //          PARAMPACK::CompareSimilarity(curCurve, curTree->GetTree()[blockInd][curveInd], headInd);
  //          if (headInd == -1) continue;
  //          /*add parent curve to child curve and delete*/
  //          PARAMPACK::ConcatTracedCurves(curCurve, curTree->GetTree()[blockInd][curveInd], headInd);
  //          curTree->GetTree()[blockInd][curveInd].clear();
  //          blockIndList.push_back(blockInd);
  //      }
  //      //clear empty curves
  //      for (auto blockIndIt : blockIndList) {
  //          auto begIt = curTree->GetTree()[blockIndIt].begin();
  //          //auto endIt = curTree->GetTree()[blockIndIt].end();
  //          for (auto it = begIt; it != curTree->GetTree()[blockIndIt].end();) {
  //              if (it->empty()) it = curTree->GetTree()[blockIndIt].erase(it);
  //              else ++it;
  //          }
  //      }
  //      //clear empty block curves
  //      //cuz vector<vector> the element is point, so erase operation cost little
  //      auto it = curTree->GetTree().begin();
  //      while (it != curTree->GetTree().end()) {
  //          if ((*it).empty()) it = curTree->GetTree().erase(it);
  //          else ++it;
  //      }
  //  }

  //  void FindNearestCurvesEnd(IDataPointer m_Source, const VectorVec5d &curCurve, int &blockInd, int &curveInd)
  //  {
  //      blockInd = curveInd = -1;
  //      NG_CREATE_DYNAMIC_CAST(NeuronPopulation, curTree, m_Source);
  //      double minDist = 100000.0, tmpDist;
  //      for (size_t i = 0; i < curTree->GetTree().size(); ++i) {
  //          for (size_t j = 0; j < curTree->GetTree()[i].size(); ++j) {
  //              if (&curCurve == &curTree->GetTree()[i][j] || curTree->GetTree()[i][j].empty()) continue;//self
  //              /*compare head of curCurve and tail of other curves*/
  //              const Vec5d &curNode = curCurve.front();
  //              tmpDist = (curNode.block(0, 0, 3, 1) - curTree->GetTree()[i][j].back().block(0, 0, 3, 1)).norm();
  //              if (tmpDist < minDist && tmpDist < 20.0) {
  //                  minDist = tmpDist;
  //                  blockInd = int(i);
  //                  curveInd = int(j);
  //              }
  //          }
  //      }
  //  }

    void CompareSimilarity(const VectorVec5d &curCurve, const VectorVec5d &dstCurve, size_t &headInd)
    {
        /*compare head of curCurve and tail of other curves*/
        const Vec5d &dstNode = dstCurve.back();
        double dist;
        headInd = std::numeric_limits<size_t>::max();
        size_t dstSz = dstCurve.size() - 1lu;// , curSz = curCurve.size() - 1lu;
        size_t maxLen = std::min(curCurve.size() / 3lu, dstSz);
        for (size_t i = 0; i < maxLen; ++i) {
            dist = (curCurve[i].block(0, 0, 3, 1) - dstNode.block(0, 0, 3, 1)).norm();
            if (dist < 3.0) {
                double avgDist = 0.0;
                for (size_t j = 0; j <= i; ++j) avgDist += (curCurve[i - j].block(0, 0, 3, 1) - dstCurve[dstSz - j].block(0, 0, 3, 1)).norm();
                avgDist /= double(i + 1);
                if (avgDist < 3.0) {
                    headInd = i;
                    break;
                }
            }
        }
    }

    void ConcatTracedCurves(VectorVec5d &childCurve, VectorVec5d &concatCurve, int headInd)
    {
        concatCurve.reserve(concatCurve.size() + childCurve.size() - headInd);
        std::copy(childCurve.begin() + headInd, childCurve.end(), std::back_inserter(concatCurve));
        concatCurve.swap(childCurve);
    }

}
