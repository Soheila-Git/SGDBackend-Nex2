import { fromJS } from 'immutable';

// temp fixture
const DEFAULT_STATE = fromJS({
  activeEntries: [
    {
      id: '12345abc',
      title: 'Lorem Ipsum dalor it Clylin Dependent Protein Serine',
      author: 'Lorem et al.',
      journal: 'Nucleic Acids Research',
      abstract: 'Sed ut perspiciatis unde omnis iste natus error sit voluptatem accusantium doloremque laudantium, totam rem aperiam, eaque ipsa quae ab illo inventore veritatis et quasi architecto beatae vitae dicta sunt explicabo.',
      status: 'reviewing'
    },
    {
      id: '67990cde',
      title: 'Lorem Ipsum dalor it Clylin Dependent Protein Serine',
      author: 'Lorem et al.',
      journal: 'Nucleic Acids Research',
      abstract: 'Sed ut perspiciatis unde omnis iste natus error sit voluptatem accusantium doloremque laudantium, totam rem aperiam, eaque ipsa quae ab illo inventore veritatis et quasi architecto beatae vitae dicta sunt explicabo.',
      status: 'reviewing'
    },
    {
      id: '67990cde',
      title: 'Lorem Ipsum dalor it Clylin Dependent Protein Serine',
      author: 'Lorem et al.',
      journal: 'Nucleic Acids Research',
      abstract: 'Sed ut perspiciatis unde omnis iste natus error sit voluptatem accusantium doloremque laudantium, totam rem aperiam, eaque ipsa quae ab illo inventore veritatis et quasi architecto beatae vitae dicta sunt explicabo.',
      status: 'reviewing'
    }
  ],
  activeLitEntry: {
    id: '12345abc',
    title: 'Yeast RAD2, a homolog of human XPG, plays a key role in the regulation of the cell cycle and actin dynamics. Biol Open',
    author: 'Lorem et al.',
    citation: 'Kang MS, et al. (2013) Yeast RAD2, a homolog of human XPG, plays a key role in the regulation of the cell cycle and actin dynamics. Biol Open',
    journal: 'Nucleic Acids Research',
    abstract: 'Sed ut perspiciatis unde omnis iste natus error sit voluptatem accusantium doloremque laudantium, totam rem aperiam, eaque ipsa quae ab illo inventore veritatis et quasi architecto beatae vitae dicta sunt explicabo.',
    status: 'reviewing'
  }
});

export default function litReducer(state = DEFAULT_STATE, action) {
  switch (action.type) {
  default:
    return state;
  }
}
