import _ from 'underscore';

const DEFAULT_STATE = {
  isAuthenticated: true,
  isAuthenticating: false,
  email: 'user123',
  csrfToken: null,
  loginError: null // can be { message: 'user not found ', error: new Errow('user not found') }
};

export default function authReducer(_state, action) {
  if (typeof _state === 'undefined') return DEFAULT_STATE;
  let state = _.clone(_state);
  switch (action.type) {
    case 'SET_CSRF_TOKEN':
      state.csrfToken = action.payload;
      return state;
      break;
    case 'RECEIVE_AUTH_RESPONSE':
      state.isAuthenticated = true;
      state.isAuthenticating = false;
      state.email = action.payload.email;
      return state;
      break;
    case 'LOGOUT':
      state.isAuthenticated = false;
      state.isAuthenticating = false;
      state.email = null;
      return state;
      break;
    default:
      return state;
  }
};
