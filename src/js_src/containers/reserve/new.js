import React, { Component } from 'react';

import FlexiForm from '../../components/forms/flexiForm';
import t from 'tcomb-form';

const TARGET_URL = '/reserve';

class GeneNameReservation extends Component {
  constructor(props) {
    super(props);
    this.state = {
      isSuccess: false
    };
  }

  renderSuccess() {
    return (
      <div>
        <p>Thanks for submitting your gene name reservation! SGD curators will review and be in touch.</p>
      </div>
    );
  }

  render() {
    if (this.state.isSuccess) {
      return this.renderSuccess();
    }
    let Author = t.struct({
      first_name: t.maybe(t.String),
      last_name: t.maybe(t.String),
      orcid: t.maybe(t.String)
    });
    let reserveSchema = t.struct({
      new_gene_name: t.maybe(t.String),
      orf_name: t.maybe(t.String),
      description: t.maybe(t.String),
      notes: t.maybe(t.String),
      first_name: t.maybe(t.String),
      last_name: t.maybe(t.String),
      email: t.maybe(t.String),
      phone_number: t.maybe(t.String),
      position: t.maybe(t.String),
      institution: t.maybe(t.String),
      profession: t.maybe(t.String),
      publication_title: t.maybe(t.String),
      journal: t.maybe(t.String),
      year: t.maybe(t.String),
      authors: t.maybe(t.list(Author))
    });
    let authorLayout = locals => {
      return (
        <div className='row'>
          <div className='column small-4'>{locals.inputs.first_name}</div>
          <div className='column small-4'>{locals.inputs.last_name}</div>
          <div className='column small-2'>{locals.inputs.orcid}</div>
          <div className='column small-2'>{locals.inputs.removeItem}</div>
        </div>
      );
    };
    let formLayout = locals => {
      return (
        <div>
          <p>* indicates required field</p>
          <p><b>Gene Name Information</b></p>
          <div className='row'>
            <div className='column small-6'>{locals.inputs.new_gene_name}</div>
            <div className='column small-6'>{locals.inputs.orf_name}</div>
          </div>
          <div>{locals.inputs.description}</div>
          <div>{locals.inputs.notes}</div>
          <p><b>Your Information</b></p>
          <div className='row'>
            <div className='column small-3'>{locals.inputs.first_name}</div>
            <div className='column small-3'>{locals.inputs.last_name}</div>
            <div className='column small-3'>{locals.inputs.email}</div>
            <div className='column small-3'>{locals.inputs.phone_number}</div>
          </div>
          <div className='row'>
            <div className='column small-4'>{locals.inputs.position}</div>
            <div className='column small-4'>{locals.inputs.institution}</div>
            <div className='column small-4'>{locals.inputs.profession}</div>
          </div>
          <p><b>Publication Information</b></p>
          <div className='row'>
            <div className='column small-6'>{locals.inputs.publication_title}</div>
            <div className='column small-4'>{locals.inputs.journal}</div>
            <div className='column small-2'>{locals.inputs.year}</div>
          </div>
          <span><a href='https://orcid.org/register' target='_new'><i className='fa fa-question-circle' /> Register for an ORCID iD</a></span>
          <div>{locals.inputs.authors}</div>
        </div>
      );
    };
    let reserveOptions = {
      template: formLayout,
      fields: {
        new_gene_name: {
          label: 'Proposed Gene Name *'
        },
        description: {
          label: 'Description of Gene Name Acronym *'
        },
        orf_name: {
          label: 'ORF Name'
        },
        first_name: {
          label: 'First Name *'
        },
        last_name: {
          label: 'Last Name *'
        },
        email: {
          label: 'Email *'
        },
        year: {
          label: 'Year *'
        },
        authors: {
          disableOrder: true,
          disableRemove: true,
          item: {
            template: authorLayout
          }
        }
      }
    };
    let _onSuccess = () => {
      this.setState({ isSuccess: true });
    };
    let _defaultData = { authors: [{ first_name: ''}] };
    t.form.Form.i18n = {
      optional: '',
      required: '',
      add: 'Add another author',
      remove: 'Remove this author'
    };
    return (
      <div>
        <h1>Reserve a Gene Name</h1>
        <FlexiForm defaultData={_defaultData} tFormOptions={reserveOptions} tFormSchema={reserveSchema} onSuccess={_onSuccess} requestMethod='POST' submitText='Send gene name reservation' updateUrl={TARGET_URL} />
      </div>
    );
  }
}

export default GeneNameReservation;