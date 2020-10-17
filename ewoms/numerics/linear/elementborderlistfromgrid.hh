// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the eWoms project.

  eWoms is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  eWoms is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with eWoms.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 * \copydoc Ewoms::Linear::ElementBorderListFromGrid
 */
#ifndef EWOMS_ELEMENT_BORDER_LIST_FROM_GRID_HH
#define EWOMS_ELEMENT_BORDER_LIST_FROM_GRID_HH

#include "overlaptypes.hh"
#include "blacklist.hh"

#include <ewoms/common/unused.hh>

#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>
#include <dune/common/version.hh>

#include <algorithm>

namespace Ewoms {
namespace Linear {
/*!
 * \brief Uses communication on the grid to find the initial seed list
 *        of indices for methods which use element-based degrees of
 *        freedom.
 */
template <class GridView, class ElementMapper>
class ElementBorderListFromGrid
{
    using PeerBlackListedEntry = BlackList::PeerBlackListedEntry;
    using PeerBlackList = BlackList::PeerBlackList;
    using PeerBlackLists = BlackList::PeerBlackLists;

    using Element = typename GridView::template Codim<0>::Entity;

    class BorderListHandle_
        : public Dune::CommDataHandleIF<BorderListHandle_, ProcessRank>
    {
    public:
        BorderListHandle_(const GridView& gridView,
                          const ElementMapper& map,
                          BlackList& blackList,
                          BorderList& borderList)
            : gridView_(gridView)
            , map_(map)
            , blackList_(blackList)
            , borderList_(borderList)
        {
            auto elemIt = gridView_.template begin<0>();
            const auto& elemEndIt = gridView_.template end<0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                if (elemIt->partitionType() != Dune::InteriorEntity) {
                    Index elemIdx = static_cast<Index>(map_.index(*elemIt));
                    blackList_.addIndex(elemIdx);
                }
            }
        }

        // data handle methods
        bool contains(int dim EWOMS_UNUSED, int codim) const
        { return codim == 0; }

#if DUNE_VERSION_NEWER(DUNE_GRID, 2,7)
        bool fixedSize(int dim EWOMS_UNUSED, int codim EWOMS_UNUSED) const
        { return true; }
#else
        bool fixedsize(int dim EWOMS_UNUSED, int codim EWOMS_UNUSED) const
        { return true; }
#endif

        template <class EntityType>
        size_t size(const EntityType& e EWOMS_UNUSED) const
        { return 2; }

        template <class MessageBufferImp, class EntityType>
        void gather(MessageBufferImp& buff, const EntityType& e) const
        {
            unsigned myIdx = static_cast<unsigned>(map_.index(e));
            buff.write(static_cast<unsigned>(gridView_.comm().rank()));
            buff.write(myIdx);
        }

        template <class MessageBufferImp>
        void scatter(MessageBufferImp& buff,
                     const Element& e,
                     size_t n EWOMS_UNUSED)
        {
            // discard the index if it is not on the process boundary
            bool isInteriorNeighbor = false;
            auto isIt = gridView_.ibegin(e);
            const auto& isEndIt = gridView_.iend(e);
            for (; isIt != isEndIt; ++isIt) {
                if (!isIt->neighbor())
                    continue;
                else if (isIt->outside().partitionType() == Dune::InteriorEntity) {
                    isInteriorNeighbor = true;
                    break;
                }
            }
            if (!isInteriorNeighbor)
                return;

            BorderIndex bIdx;

            bIdx.localIdx = static_cast<Index>(map_.index(e));
            {
                ProcessRank tmp;
                buff.read(tmp);
                bIdx.peerRank = tmp;
                peerSet_.insert(tmp);
            }
            {
                unsigned tmp;
                buff.read(tmp);
                bIdx.peerIdx = static_cast<Index>(tmp);
            }
            bIdx.borderDistance = 1;

            borderList_.push_back(bIdx);
        }

        // this template method is needed because the above one only works for codim-0
        // entities (i.e., elements) but the dune grid uses some code which causes the
        // compiler to invoke the scatter method for every codim...
        template <class MessageBufferImp, class EntityType>
        void scatter(MessageBufferImp& buff EWOMS_UNUSED,
                     const EntityType& e EWOMS_UNUSED,
                     size_t n EWOMS_UNUSED)
        { }

        const std::set<ProcessRank>& peerSet() const
        { return peerSet_; }

    private:
        GridView gridView_;
        const ElementMapper& map_;
        std::set<ProcessRank> peerSet_;
        BlackList& blackList_;
        BorderList& borderList_;
    };

    class PeerBlackListHandle_
        : public Dune::CommDataHandleIF<PeerBlackListHandle_, int>
    {
    public:
        PeerBlackListHandle_(const GridView& gridView,
                             const ElementMapper& map,
                             PeerBlackLists& peerBlackLists)
            : gridView_(gridView)
            , map_(map)
            , peerBlackLists_(peerBlackLists)
        {}

        // data handle methods
        bool contains(int dim EWOMS_UNUSED, int codim) const
        { return codim == 0; }

#if DUNE_VERSION_NEWER(DUNE_GRID, 2,7)
        bool fixedSize(int dim EWOMS_UNUSED, int codim EWOMS_UNUSED) const
        { return true; }
#else
        bool fixedsize(int dim EWOMS_UNUSED, int codim EWOMS_UNUSED) const
        { return true; }
#endif

        template <class EntityType>
        size_t size(const EntityType& e EWOMS_UNUSED) const
        { return 2; }

        template <class MessageBufferImp, class EntityType>
        void gather(MessageBufferImp& buff, const EntityType& e) const
        {
            buff.write(static_cast<int>(gridView_.comm().rank()));
            buff.write(static_cast<int>(map_.index(e)));
        }

        template <class MessageBufferImp, class EntityType>
        void scatter(MessageBufferImp& buff, const EntityType& e, size_t n EWOMS_UNUSED)
        {
            int peerRank;
            int peerIdx;
            Index localIdx;

            buff.read(peerRank);
            buff.read(peerIdx);
            localIdx = static_cast<Index>(map_.index(e));

            PeerBlackListedEntry pIdx;
            pIdx.nativeIndexOfPeer = static_cast<Index>(peerIdx);
            pIdx.myOwnNativeIndex = static_cast<Index>(localIdx);

            peerBlackLists_[static_cast<ProcessRank>(peerRank)].push_back(pIdx);
        }

        const PeerBlackList& peerBlackList(ProcessRank peerRank) const
        { return peerBlackLists_.at(peerRank); }

    private:
        GridView gridView_;
        const ElementMapper& map_;
        PeerBlackLists peerBlackLists_;
    };

public:
    ElementBorderListFromGrid(const GridView& gridView, const ElementMapper& map)
        : gridView_(gridView)
        , map_(map)
    {
        BorderListHandle_ blh(gridView, map, blackList_, borderList_);
        gridView.communicate(blh,
                             Dune::InteriorBorder_All_Interface,
                             Dune::BackwardCommunication);

        PeerBlackListHandle_ pblh(gridView, map, peerBlackLists_);
        gridView.communicate(pblh,
                             Dune::InteriorBorder_All_Interface,
                             Dune::BackwardCommunication);

        auto peerIt = blh.peerSet().begin();
        const auto& peerEndIt = blh.peerSet().end();
        for (; peerIt != peerEndIt; ++peerIt) {
            blackList_.setPeerList(*peerIt, pblh.peerBlackList(*peerIt));
        }
    }

    // Access to the border list.
    const BorderList& borderList() const
    { return borderList_; }

    // Access to the black-list indices.
    const BlackList& blackList() const
    { return blackList_; }

private:
    const GridView gridView_;
    const ElementMapper& map_;

    BorderList borderList_;

    BlackList blackList_;
    PeerBlackLists peerBlackLists_;
};

} // namespace Linear
} // namespace Ewoms

#endif
