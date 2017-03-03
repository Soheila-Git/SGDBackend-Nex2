import React, { Component } from 'react';
import { connect } from 'react-redux';
import Dropdown, { DropdownTrigger, DropdownContent } from 'react-simple-dropdown';

import style from './style.css';
// import fetchData from '../../lib/fetchData';
import { removeEntry } from './triageActions';
import { setMessage } from '../../actions/metaActions';
import TagList from './tagList';

// const TRIAGE_URL = '/reference/triage';
// const PROMOTE_URL_SUFFIX = 'promote';

class TriageControls extends Component {
  handleDiscardEntry(e) {
    e.preventDefault();
    let id = this.props.id;
    this.props.dispatch(removeEntry(id));
    let message = `${this.props.citation} discarded.`;
    this.props.dispatch(setMessage(message));
    // let url = `${TRIAGE_URL}/${id}`;
    // let fetchOptions = {
    //   type: 'DELETE',
    //   headers: {
    //     'X-CSRF-Token': window.CSRF_TOKEN,        
    //   }
    // };
    // fetchData(url, fetchOptions).then( () => {
    //   this.props.dispatch(promoteEntry(id));
    // });
  }

  handlePromoteEntry(e) {
    e.preventDefault();
    let id = this.props.id;
    this.props.dispatch(removeEntry(id));
    let message = `${this.props.citation} added to database.`;
    this.props.dispatch(setMessage(message));
    // let url = `${TRIAGE_URL}/${id}/${PROMOTE_URL_SUFFIX}`;
    // let fetchOptions = {
    //   type: 'PUT',
    //   headers: {
    //     'X-CSRF-Token': window.CSRF_TOKEN,        
    //   }
    // };
    // fetchData(url, fetchOptions).then( () => {
    //   this.props.dispatch(promoteEntry(id));
    // });
  }

  renderAssign() {
    let assignee = this.props.data.assignee;
    let node;
    if (assignee) {
      node = (
        <div>
          <span>{assignee}</span>
        </div>
      ); 
    } else {
      // let onClaim 
      node = (
        <div>
          <span>{assignee}</span>
        </div>
      );
    }
    return (
      <div>
        <label>Assignee</label>
        {node}
      </div>
    );
  }

  renderTags() {
    return (
      <Dropdown>
        <DropdownTrigger className='button small'>Tags <i className='fa fa-caret-down' /></DropdownTrigger>
        <DropdownContent className={`dropdownContent ${style.tagList}`}>
          <TagList id={this.props.id} />
        </DropdownContent>
      </Dropdown>
    );
  }

  render() {
    return (
      <div className='row'>
        {this.renderAssign()}
        <div className='column small-6'>
          {this.renderTags()}
        </div>
        <div className='column small-6 text-right'>
          <a className='button small' onClick={this.handlePromoteEntry.bind(this)}><i className='fa fa-check' /> Add to Database</a>
          <a className='button secondary small' onClick={this.handleDiscardEntry.bind(this)}><i className='fa fa-trash' /> Discard</a>
        </div>
      </div>
    );
  }
}

TriageControls.propTypes = {
  citation: React.PropTypes.string,
  data: React.PropTypes.object,
  dispatch: React.PropTypes.func,
  id: React.PropTypes.number
};

function mapStateToProps() {
  return {};
}

export default connect(mapStateToProps)(TriageControls);
